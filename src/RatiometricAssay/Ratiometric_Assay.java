/*
 * Copyright (C) 2017 Dave Barry <david.barry at crick.ac.uk>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package RatiometricAssay;

import MetaData.ParamsReader;
import Profile.PeakFinder;
import UtilClasses.GenUtils;
import UtilClasses.GenVariables;
import UtilClasses.Utilities;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.plugin.LutLoader;
import ij.plugin.PlugIn;
import ij.process.Blitter;
import ij.process.FloatBlitter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.LUT;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;

public class Ratiometric_Assay implements PlugIn {

    private final String title = "Rationmetric Assay Analyser", um = String.format("%c%c", IJ.micronSymbol, 'm');
    private static double spatialRes, timeRes;

    public static void main(String args[]) {
        Ratiometric_Assay instance = new Ratiometric_Assay();
        instance.run(null);
        System.exit(0);
    }

    public void run(String arg) {
        ImageStack stack1, stack2;
        if (IJ.getInstance() == null) {
            stack1 = IJ.openImage().getImageStack();
            stack2 = IJ.openImage().getImageStack();
        } else {
            ImagePlus[] imps = GenUtils.specifyInputs(new String[]{"Stack 1", "Stack 2"});
            readParamsFromImage(imps[0]);
            stack1 = imps[0].getImageStack();
            stack2 = imps[1].getImageStack();
        }
        if (stack1.size() != stack2.size()
                || stack1.getWidth() != stack2.getWidth()
                || stack1.getHeight() != stack2.getHeight()) {
            GenUtils.error("Stack dimensions must match.");
            return;
        }
        if (!showDialog()) {
            return;
        }
        String resultsDir;
        try {
            resultsDir = GenUtils.openResultsDirectory(String.format("%s%s%s",
                    (Utilities.getFolder(new File(IJ.getDirectory("current")),
                            "Specify directory for output files...",
                            true)).getAbsolutePath(), File.separator, title));
        } catch (Exception e) {
            GenUtils.error(e.toString());
            return;
        }
        ImageStack output = new ImageStack(stack1.getWidth(), stack1.getHeight());
        int n = stack1.getSize();
        for (int i = 1; i <= n; i++) {
            FloatProcessor slice1 = stack1.getProcessor(i).convertToFloatProcessor();
            FloatProcessor slice2 = stack2.getProcessor(i).convertToFloatProcessor();
            FloatProcessor outSlice = (FloatProcessor) slice1.duplicate();
            (new FloatBlitter(outSlice)).copyBits(slice2, 0, 0, Blitter.DIVIDE);
            output.addSlice(outSlice);
        }
        ImagePlus outputImp = new ImagePlus("Ratiometric Output", output);
        outputImp.resetDisplayRange();
        LUT sixteenColors;
        if (IJ.getInstance() != null) {
            sixteenColors = LutLoader.openLut(String.format("%s%s%s%s", IJ.getDirectory("imagej"), "luts", File.separator, "16_colors.lut"));
        } else {
            sixteenColors = LutLoader.openLut("C:\\Users\\barryd\\fiji-nojre\\Fiji.app\\luts\\16_colors.lut");
        }
        outputImp.setLut(sixteenColors);
        outputImp.show();
        IJ.saveAs(outputImp, "TIF", String.format("%s%s%s", resultsDir, File.separator, "output.tif"));
        try {
            double[][] data = compressSlices(output, resultsDir);
            plotProfilePoints(data, resultsDir);
        } catch (Exception e) {
            GenUtils.error(e.toString());
            return;
        }
        IJ.showStatus(String.format("%s done.", title));
    }

    double[][] compressSlices(ImageStack stack, String dir) throws IOException, FileNotFoundException {
        File output = new File(String.format("%s%s%s", dir, File.separator, "Map.csv"));
        CSVPrinter printer = new CSVPrinter(new OutputStreamWriter(new FileOutputStream(output), GenVariables.ISO), CSVFormat.EXCEL);
        int n = stack.size();
        int width = stack.getWidth();
        int height = stack.getHeight();
        double[][] data = new double[n][width];
        for (int i = 1; i <= n; i++) {
            ImageProcessor slice = stack.getProcessor(i);
            for (int x = 0; x < width; x++) {
                double sum = 0.0;
                for (int y = 0; y < height; y++) {
                    sum += slice.getPixelValue(x, y);
                }
                data[i - 1][x] = sum / height;
                printer.print(String.valueOf(data[i - 1][x]));
            }
            printer.println();
        }
        printer.close();
        return data;
    }

    void plotProfilePoints(double[][] data, String dir) throws IOException, FileNotFoundException {
        File output = new File(String.format("%s%s%s", dir, File.separator, "ProfilePoints.csv"));
        CSVPrinter printer = new CSVPrinter(new OutputStreamWriter(new FileOutputStream(output), GenVariables.ISO), CSVFormat.EXCEL);
        printer.printRecord("Time (s)", "Scratch Position " + um, "Peak Position" + um, "Edge Position " + um);
        for (int i = 0; i < data.length; i++) {
            printer.print(i * timeRes);
            int[] indices = PeakFinder.findMaxAndSides(PeakFinder.smoothData(data[i], 10.0));
            for (int j = 0; j < indices.length; j++) {
                if (indices[j] > 0) {
                    printer.print(indices[j] * spatialRes);
                } else {
                    printer.print("not found");
                }
            }
            printer.println();
        }
        printer.close();
    }

    boolean showDialog() {
        GenericDialog gd = new GenericDialog(title);
        gd.addNumericField("Spatial Resolution:", spatialRes, 3, 5, um);
        gd.addNumericField("Frame Interval:", timeRes, 3, 5, "seconds");
        gd.showDialog();
        if (gd.wasCanceled()) {
            return false;
        }
        spatialRes = gd.getNextNumber();
        timeRes = gd.getNextNumber();
        if (gd.wasOKed()) {
            return true;
        } else {
            return false;
        }
    }

    void readParamsFromImage(ImagePlus imp) {
        ParamsReader reader = new ParamsReader(imp);
        spatialRes = reader.getSpatialRes();
        timeRes = reader.getFrameRate();
    }
}
