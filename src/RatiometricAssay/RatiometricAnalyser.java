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

import Binary.BinaryMaker;
import Binary.EDMMaker;
import ImageProcessing.ImageBlurrer;
import MetaData.ParamsReader;
import Profile.PeakFinder;
import UtilClasses.GenUtils;
import UtilClasses.GenVariables;
import UtilClasses.Utilities;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.plugin.LutLoader;
import ij.process.AutoThresholder;
import ij.process.Blitter;
import ij.process.ByteBlitter;
import ij.process.ByteProcessor;
import ij.process.FloatBlitter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.LUT;
import ij.process.StackStatistics;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

public class RatiometricAnalyser {

    private final String title = "Ratiometric Assay Analyser", um = String.format("%c%c", IJ.micronSymbol, 'm');
    private static double spatialRes, timeRes, threshold = 0.5;
    public static final String[] THRESH_METHODS = AutoThresholder.getMethods();
    private static String threshMethod = AutoThresholder.Method.Triangle.toString();
    private static double maskBlurRadius = 1.0;
    private static double sigBlurRadius = 20.0;
    private static int holeSize = 10;

    public void analyse(ImageStack stack1, ImageStack stack2) {
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
        ImageStack edmStack = makeEDMStack(stack1.duplicate(), resultsDir);
        plotVelocities(estimateVelocity(edmStack));
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
            sixteenColors = LutLoader.openLut("C:\\Users\\barryd\\FIJI\\DEV\\Fiji.app\\luts\\16_colors.lut");
        }
        outputImp.setLut(sixteenColors);
        outputImp.show();
        IJ.saveAs(outputImp, "TIF", String.format("%s%s%s", resultsDir, File.separator, "output.tif"));
        try {
            double[][] data = indexWithDistanceMap(output, edmStack, resultsDir);
            plotProfilePoints(data, resultsDir);
        } catch (Exception e) {
            GenUtils.error(e.toString());
            return;
        }
        IJ.showStatus(String.format("%s done.", title));
    }

    public static ImageStack makeBinaryStack(ImageStack input, double maskBlurRadius, String threshMethod, int holeSize) {
        ImageBlurrer.blurStack(input, maskBlurRadius);
        return BinaryMaker.makeBinaryStack(new ImagePlus("", input), threshMethod, -1, holeSize);
    }

    ImageStack makeEDMStack(ImageStack input, String directory) {
        ImageStack binary = makeBinaryStack(input, maskBlurRadius, threshMethod, holeSize);
        ImageStack output = EDMMaker.makeEDM(binary);
        IJ.saveAs(new ImagePlus("", output), "TIF", String.format("%s%s%s", directory, File.separator, "EDM"));
        return output;
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

    double[][] indexWithDistanceMap(ImageStack stack, ImageStack edmStack, String dir) throws IOException, FileNotFoundException {
        File output = new File(String.format("%s%s%s", dir, File.separator, "Map.csv"));
        CSVPrinter printer = new CSVPrinter(new OutputStreamWriter(new FileOutputStream(output), GenVariables.ISO), CSVFormat.EXCEL);
        int n = stack.size();
        StackStatistics stats = new StackStatistics(new ImagePlus("", edmStack));
        int maxDist = (int) Math.round(stats.max);
        int width = stack.getWidth();
        int height = stack.getHeight();
        int imArea = width * height;
        double[][] data = new double[n][maxDist + 1];
        for (int i = 1; i <= n; i++) {
            ImageProcessor slice = stack.getProcessor(i).convertToFloatProcessor();
            slice.multiply(1.0 / imArea);
            ImageProcessor edmSlice = edmStack.getProcessor(i);
            for (int x = 0; x < width; x++) {
                for (int y = 0; y < height; y++) {
                    int edmDist = (int) Math.round(edmSlice.getPixelValue(x, y));
                    data[i - 1][edmDist] = slice.getPixelValue(x, y);
                }
            }
            for (int j = 0; j < maxDist; j++) {
                printer.print(String.valueOf(data[i - 1][j]));
            }
            printer.println();
        }
        printer.close();
        return data;
    }

    void plotProfilePoints(double[][] data, String dir) throws IOException, FileNotFoundException {
        File output = new File(String.format("%s%s%s", dir, File.separator, "ProfilePoints.csv"));
        CSVPrinter printer = new CSVPrinter(new OutputStreamWriter(new FileOutputStream(output), GenVariables.ISO), CSVFormat.EXCEL);
        printer.printRecord("Time (s)", "x1 " + um, "x2 " + um, "Width " + um);
        double[][] smoothedData = PeakFinder.smoothData(data, sigBlurRadius);
        double[] extrema = PeakFinder.findMinAndMax(smoothedData);
        for (int i = 0; i < data.length; i++) {
            printer.print(i * timeRes);
            int[] indices = PeakFinder.findRegionWidth(smoothedData[i], threshold, extrema[0], extrema[1]);
            for (int j = 0; j < indices.length; j++) {
                if (indices[j] > 0) {
                    printer.print(indices[j] * spatialRes);
                } else {
                    printer.print("not found");
                }
            }
            if (indices[0] > 0 && indices[1] > 0) {
                printer.print((indices[1] - indices[0]) * spatialRes);
            }
            printer.println();
        }
        printer.close();
    }

    boolean showDialog() {
        GenericDialog gd = new GenericDialog(title);
        gd.addNumericField("Spatial Resolution:", spatialRes, 3, 5, um);
        gd.addNumericField("Frame Interval:", timeRes, 3, 5, "seconds");
        gd.addNumericField("Active Threshold:", threshold, 3, 5, "");
        gd.showDialog();
        if (gd.wasCanceled()) {
            return false;
        }
        spatialRes = gd.getNextNumber();
        timeRes = gd.getNextNumber();
        threshold = gd.getNextNumber();
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

    private double[] estimateVelocity(ImageStack edm) {
        ImageStack binaryStack = BinaryMaker.makeBinaryStack(new ImagePlus("", edm.duplicate()), null, 1, holeSize);
        int size = binaryStack.size();
        double[] vels = new double[size];
        vels[0] = 0.0;
        SummaryStatistics stats = new SummaryStatistics();
        for (int i = 2; i <= size - 1; i++) {
            ByteProcessor slice1 = (ByteProcessor) binaryStack.getProcessor(i - 1).duplicate();
            ByteProcessor slice2 = (ByteProcessor) binaryStack.getProcessor(i + 1).duplicate();
            ByteBlitter blitter = new ByteBlitter(slice2);
            blitter.copyBits(slice1, 0, 0, Blitter.DIFFERENCE);
            double area = slice2.getStatistics().histogram[255];
            ImageProcessor edmslice = edm.getProcessor(i);
            double perim = edmslice.getStatistics().histogram[1];
            vels[i - 1] = (spatialRes * area) / (2.0 * perim * timeRes); // Divide by 2 because velocties are calculated over two frame intervals
            stats.addValue(vels[i - 1]);
        }
        vels[size - 1] = vels[size - 2];
        IJ.log(String.format("Estimated mean cell velocity: %f %s", stats.getMean(), um));
        return vels;
    }

    private void plotVelocities(double[] velocities) {
        double[] timePoints = new double[velocities.length];
        for (int i = 0; i < velocities.length; i++) {
            timePoints[i] = i * timeRes;
        }
        Plot plot = new Plot("Instantaneous Velocities", "Time", "Instantaneous Velocity", timePoints, velocities);
        plot.show();
    }
}
