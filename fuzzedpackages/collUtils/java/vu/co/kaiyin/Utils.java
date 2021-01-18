package vu.co.kaiyin;

import org.apache.commons.io.FilenameUtils;

import java.io.*;

/**
 * Created by IDEA on 14/06/15.
 */
public class Utils {
    public static int[] seq(int start, int end) {
        if(end < start) {
            end = start;
        }
        int[] res = new int[end - start + 1];
        for(int i = start; i <= end; i++) {
            res[i - start] = i;
        }
        return res;
    }

    public static int countLines(String filename) throws IOException {
        BufferedReader reader = new BufferedReader(
                new FileReader(filename)
        );
        int lines = 0;
        try{
            while(reader.readLine() != null) {
                lines++;
            }
        } finally {
            reader.close();
        }
        return lines;
    }

    public static String shiftedFilename(String filename, int nShift) {
        String ext = FilenameUtils.getExtension(filename);
        return String.format("%s_shift_%04d.%s", FilenameUtils.removeExtension(filename), nShift, ext);
    }

    public static void printMat(int[][] x, int startRow, int endRow, int startColumn, int endColumn) {
        if(x.length == 0 || x[0].length == 0) return;
        if(startRow < 0) startRow = 0;
        if(startRow >= x.length) startRow = x.length - 1;
        if(endRow < startRow) endRow = startRow + 1;
        if(endRow > x.length) endRow = x.length;
        if(startColumn < 0) startColumn = 0;
        if(startColumn >= x[0].length) startColumn = x[0].length - 1;
        if(endColumn < startColumn) endColumn = startColumn + 1;
        if(endColumn > x[0].length) endColumn = x[0].length;
        for(int i = startRow; i < endRow; i++) {
            for(int j = startColumn; j < endColumn; j++) {
                System.out.printf("%4d", x[i][j]);
            }
            System.out.println();
        }
        System.out.println();
    }

    public static void truncateFromEnd(File filename, int n) throws Exception {
        if(n < 0) {
            throw new Exception("Can't truncate by a negative number");
        }
        try (RandomAccessFile raf = new RandomAccessFile(filename, "rw")) {
            long originalLength = raf.length();
            long newLength = originalLength - (long) n;
            if(newLength < 0) {
                newLength = 0;
            }
            raf.setLength(newLength);
        }
    }

    public static void truncateFromEnd(String filename, int n) throws Exception {
        truncateFromEnd(new File(filename), n);
    }

    public static void printMat(int[][] x) {
        printMat(x, 0, x.length, 0, x[0].length);
    }

    public static void main(String[] args) throws Exception {
//        long t1, t2;
//        t1 = System.currentTimeMillis();
//        System.out.println(countLines("/tmp/big.txt"));
//        t2 = System.currentTimeMillis();
//        System.out.println(t2 - t1);
//        FileWriter fileWriter = new FileWriter("/tmp/out.txt");
//        fileWriter.write("a\nb");
//        fileWriter.close();
//        System.out.println(countLines("/tmp/out.txt"));
        FileOutputStream fo = new FileOutputStream("/tmp/bin.out");
        byte[] data = new byte[] {1, 2, 3, 4, 5};
        fo.write(data);
        fo.close();
        truncateFromEnd("/tmp/bin.out", 1);
    }
}
