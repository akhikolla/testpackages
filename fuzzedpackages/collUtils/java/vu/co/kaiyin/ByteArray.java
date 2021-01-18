package vu.co.kaiyin;

import java.util.Arrays;

public class ByteArray {
    private final byte[] array;

    public ByteArray(byte[] array) {
        this.array = array;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        return Arrays.equals(array, ((ByteArray) o).array);
    }

    @Override
    public int hashCode() {
        return Arrays.hashCode(array);
    }

    @Override
    public String toString() {
        return Arrays.toString(array);
    }
}
