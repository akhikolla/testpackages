package vu.co.kaiyin;

import java.security.InvalidParameterException;

/**
 * Created by IDEA on 13/06/15.
 */
public class ShiftCodes {
    private final byte[][] shiftMatrix;

    public ShiftCodes(int[][] collapseMatrix) {
        shiftMatrix = new byte[256][256];
        for (int i = 0; i < 256; i++) {
            for (int j = 0; j < 256; j++) {
                shiftMatrix[i][j] = GenoBytes.genoByte(GenoBytes.collapse(
                        (byte) i,
                        (byte) j,
                        collapseMatrix));
            }
        }
    }

    public ShiftCodes() {
        this(GenoBytes.defaultCollpaseMatrix);
    }

    public byte lookup(byte b1, byte b2) {
        return shiftMatrix[b1 & 0xff][b2 & 0xff];
    }

    public byte[] lookup(byte[] bs1, byte[] bs2) {
        if(bs1.length != bs2.length) {
            throw new InvalidParameterException("bs1 and bs2 must be of the same length");
        }
        byte[] res = new byte[bs1.length];
        for(int i = 0; i < bs1.length; i++) {
            res[i] = lookup(bs1[i], bs2[i]);
        }
        return res;
    }

    public static void main(String[] args) {
        ShiftCodes shiftCodes1 = new ShiftCodes();
        ShiftCodes shiftCodes2 = new ShiftCodes(GenoBytes.antidiangonal);
        System.out.println(shiftCodes1.lookup((byte) 0b11, (byte) 0b01));
    }
}
