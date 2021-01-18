package vu.co.kaiyin;

import java.io.FileOutputStream;
import java.io.IOException;
import java.security.InvalidParameterException;

/**
 * Created by IDEA on 13/06/15.
 */
public class GenoBytes {

    public static int[][] defaultCollpaseMatrix = new int[][] {
            {0, 0, 0, 0},
            {0, 1, 1, 1},
            {0, 1, 0, 3},
            {0, 1, 3, 3}
    };

    public static int[][] antidiangonal = new int[][] {
            {3, 1, 3, 0},
            {1, 1, 1, 1},
            {3, 1, 0, 3},
            {0, 1, 3, 3}
    };

    public static byte genoByte(byte b1, byte b2, byte b3, byte b4) {
        if(! (
                rightByte(b1) &&
                        rightByte(b2) &&
                        rightByte(b3) &&
                        rightByte(b4)
        )) {
            throw new InvalidParameterException("genotype should be between 0 and 3 (inclusive)");
        }
        return (byte) (
                (b1 << 6)
                        | (b2 << 4)
                        | (b3 << 2)
                        | b4);
    }

    public static byte genoByte(byte[] bytes) {
        if(bytes.length != 4) {
            throw new InvalidParameterException("genos must be of length 4");
        }
        return genoByte(bytes[0], bytes[1], bytes[2], bytes[3]);
    }

    public static byte genoByte(int b1, int b2, int b3, int b4) {
        return genoByte((byte) b1, (byte) b2, (byte) b3, (byte) b4);
    }

    public static byte genoByte(int[] genos) {
        byte[] bytes = new byte[4];
        for(int i = 0; i < 4; i++) {
            bytes[i] = (byte) genos[i];
        }
        return genoByte(bytes);
    }

    public static boolean rightByte(byte b) {
        return b >= 0 && b < 4;
    }

    public static boolean rightGeno(int g) {
        return g >= 0 && g < 4;
    }

    public static int[] byteGeno(byte b) {
        int i1 = (int) ((b & 0b11000000) >> 6);
        int i2 = (int) ((b & 0b00110000) >> 4);
        int i3 = (int) ((b & 0b00001100) >> 2);
        int i4 = (int) (b & 0b00000011);
        return new int[] {i1, i2, i3, i4};
    }

    public static int collapse(int geno1, int geno2, int[][] collpaseMatrix) {
        if(collpaseMatrix.length != 4 || collpaseMatrix[0].length != 4
                || collpaseMatrix[1].length != 4
                || collpaseMatrix[2].length != 4
                || collpaseMatrix[3].length != 4
                ) {
            throw new InvalidParameterException("collapse matrix must be 4x4");
        }
        if(! (rightGeno(geno1) && rightGeno(geno2))) {
            throw new InvalidParameterException("genotype must be between 0 and 3 (inclusive)");
        }
        return collpaseMatrix[geno1][geno2];
    }


    public static int collapse(int geno1, int geno2) {
        return collapse(geno1, geno2, defaultCollpaseMatrix);
    }

    public static int[] collapse(int[] genos1, int[] genos2, int[][] collapseMatrix) {
        if(genos1.length != genos2.length) {
            throw new InvalidParameterException("genos1 and genos2 must be of the same length");
        }
        int[] gCollapsed = new int[genos1.length];
        for(int i = 0; i < genos1.length; i++) {
            gCollapsed[i] = collapse(genos1[i], genos2[i], collapseMatrix);
        }
        return gCollapsed;
    }

    public static int[] collapse(int[] genos1, int[] genos2) {
        return collapse(genos1, genos2, defaultCollpaseMatrix);
    }

    public static int[] collapse(byte b1, byte b2) {
        return collapse(byteGeno(b1), byteGeno(b2));
    }

    public static int[] collapse(byte b1, byte b2, int[][] collapseMatrix) {
        return collapse(byteGeno(b1), byteGeno(b2), collapseMatrix);
    }

    public static byte[] collapse(byte[] bs1, byte[] bs2, int[][] collapseMatrix)  {
        if(bs1.length != bs2.length) {
            throw new InvalidParameterException("bs1 and bs2 must be of the same length");
        }
        byte[] collapsedBytes = new byte[bs1.length];
        for(int i = 0; i < bs1.length; i++) {
            collapsedBytes[i] = genoByte(collapse(bs1[i], bs2[i], collapseMatrix));
        }
        return collapsedBytes;
    }

    public static byte[] collapse(byte[] bs1, byte[] bs2) {
        return collapse(bs1, bs2, defaultCollpaseMatrix);
    }




    public static void main(String[] args) throws IOException {
//        byte b = shiftcodes.genoByte((byte) 0x01, (byte) 0x03, (byte) 0x00, (byte) 0x02);
//        byte b = shiftcodes.genoByte((byte) 0b01, (byte) 0b11, (byte) 0b00, (byte) 0b10);
        byte b = GenoBytes.genoByte(new int[] {1, 3, 0, 2});
        FileOutputStream fileOutputStream = new FileOutputStream("/tmp/x.bin");
        fileOutputStream.write(new byte[] {b});
        System.out.println(((byte) 0x11) > 0);
    }

}
