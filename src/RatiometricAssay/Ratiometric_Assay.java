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

import Profile.PeakFinder;
import UtilClasses.GenUtils;
import UtilClasses.GenVariables;
import UtilClasses.Utilities;
import ij.IJ;
import ij.ImageStack;
import ij.plugin.PlugIn;
import ij.process.Blitter;
import ij.process.FloatBlitter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;

public class Ratiometric_Assay implements PlugIn {

    private final String title = "Rationmetric Assay Analyser";

    public static void main(String args[]) {
        Ratiometric_Assay instance = new Ratiometric_Assay();
        instance.run(null);
        System.exit(0);
    }

    public void run(String arg) {
        ImageStack stack1 = IJ.openImage().getImageStack();
        ImageStack stack2 = IJ.openImage().getImageStack();
        if (stack1.size() != stack2.size()
                || stack1.getWidth() != stack2.getWidth()
                || stack1.getHeight() != stack2.getHeight()) {
            GenUtils.error("Stack dimensions must match.");
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
        try {
            double[][] data = compressSlices(output, resultsDir);
            plotProfilePoints(data, resultsDir);
        } catch (Exception e) {
            GenUtils.error(e.toString());
            return;
        }
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
        printer.printRecord("Scratch", "Peak", "Edge");
        for (int i = 0; i < data.length; i++) {
            int[] indices = PeakFinder.findMaxAndSides(data[i], 10.0);
            Integer[]  ii = new Integer[indices.length];
            for(int j=0;j<indices.length;j++){
                ii[j]=indices[j];
            }
            printer.printRecord(Arrays.asList(ii));
        }
        printer.close();
    }
}

//folder_name1 = uigetdir(pwd,'Select Folder 1');
//folder_name2 = uigetdir(pwd,'Select Folder 2');
//
//list1 = dir(folder_name1);
//list2 = dir(folder_name2);
//
//N1 = size(list1);
//N2 = size(list2);
//
//if(N1 ~= N2)
//    print('Channels must contain equal numbers of images.');
//end
//
//init = false;
//
//for f=1:N1
//    [pathstr1,name1,ext1] = fileparts(list1(f).name);
//    [pathstr2,name2,ext2] = fileparts(list2(f).name);
//    if(strcmp(ext1,'.tif') || strcmp(ext1, '.tiff'))        
//        file1 = strcat(folder_name1, filesep ,list1(f).name);
//        info1 = imfinfo(file1);
//        num_images1 = numel(info1);
//
//        file2 = strcat(folder_name2, filesep ,list2(f).name);
//        info2 = imfinfo(file2);
//        num_images2 = numel(info2);
//
//        if(num_images1 ~= num_images2)
//            print('Channels must contain equal numbers of images.');
//        end
//
//        timeRes = 300.0 / 60.0;
//        height = info1.Height;
//        width = info1.Width;
//        res = info1.XResolution;
//        X = (0:1:width-1) / res;
//        Y = (0:1:num_images1-1) * timeRes;
//
//        for i=1:num_images1
//            image1 = double(imread(file1,i));
//            image2 = double(imread(file2,i));
//            for y=1:height
//                for x=1:width
//                    ratios(i,y,x) = image1(y,x) / image2(y,x);
//                end
//            end
//        end
//        if(~init)
//            sumratios = zeros(num_images1, width);
//            init=true;
//        end
//        for i=1:num_images1
//            for y=1:height
//                for x=1:width
//                    sumratios(i,x) = sumratios(i,x) + ratios(i,y,x);
//                end
//            end
//        end
//    end
//end
