package vu.co.kaiyin;

import com.google.common.collect.BiMap;
import com.google.common.collect.HashBiMap;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;

import java.io.*;
import java.nio.file.FileAlreadyExistsException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.security.InvalidParameterException;
import java.util.Arrays;

/**
 * Created by IDEA on 16/03/15.
 */
public class Bed {
    private String bedFilename;
    private String bimFilename;
    private String famFilename;
    private int bytesSnp;
    private int nIndividuals;
    private int nIndividualsApparent;
    private long bedSize;
    private int nSNPs;
    private byte[] magicBytes;
    private BiMap<String, Integer> snpIndexMap;

    public Bed(String bedFilename, int bytesSnp, int nIndividuals) throws IOException {
        this(bedFilename);
    }

    public Bed(String bedFilename) throws IOException {
        if(! FilenameUtils.getExtension(bedFilename).equals("bed")) {
            throw new InvalidParameterException("bedFilename is not a bed file");
        }
        this.bedFilename = bedFilename;
        bimFilename = bimFilename(bedFilename);
        famFilename = famFilename(bedFilename);

        magicBytes = new byte[] {(byte) 0x6c, (byte) 0x1b, (byte) 0x01};

        nSNPs = Utils.countLines(bimFilename);
        nIndividuals = Utils.countLines(famFilename);
        bytesSnp = (int) Math.ceil(nIndividuals / 4.0);
        nIndividualsApparent = bytesSnp * 4;

        bedSize = (new File(bedFilename)).length();
        long theoBedSize = ((long) bytesSnp) * ((long) nSNPs) + ((long) 3);
        if(bedSize != theoBedSize) {
            throw new InvalidParameterException("bed file size not correct.");
        }

        snpIndexMap = HashBiMap.create();
    }

    public void readSnpIndex() throws FileNotFoundException {
        BufferedReader br = new BufferedReader(new FileReader(bimFilename));
        String line;
        String snp;
        int counter;
        try{
            for(line = br.readLine(), counter = 1; line != null; line = br.readLine(), counter++) {
                snp = line.split("\\s+")[1];
                snpIndexMap.put(snp, new Integer(counter));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public int getSnpIndex(String snp) throws FileNotFoundException {
        if(snpIndexMap.isEmpty()) {
            readSnpIndex();
        }
        return snpIndexMap.get(snp).intValue();
    }

    public int[] getSnpIndex(String[] snp) throws FileNotFoundException {
        int[] indices = new int[snp.length];
        for(int i = 0; i < snp.length; i++) {
            indices[i] = getSnpIndex(snp[i]);
        }
        return indices;
    }

    public String bimFilename(String filename) {
        return String.format("%s.bim", FilenameUtils.removeExtension(filename));
    }

    public String famFilename(String filename) {
        return String.format("%s.fam", FilenameUtils.removeExtension(filename));
    }

    public int[][] readBed(int[] snpVec, boolean transpose) throws IOException {
        int nSNPtoRead = snpVec.length;
        int numbers_per_snp = bytesSnp * 4;
        byte[] snpBuffer;
        int[][] geno = new int[nSNPtoRead][nIndividuals];
        int[][] geno_transposed = new int[nIndividuals][nSNPtoRead];
        RandomAccessFile bedFile = new RandomAccessFile(bedFilename, "r");
        try {
            // read the raw bytes
            for (int i = 0; i < nSNPtoRead; i++) {
                snpVec[i]--;
                snpBuffer = readSnpBytes(bedFile, snpVec[i]);
                // translate bytes into genotype numbers
                // every byte expanded into 4 integers:
                // http://fs2.directupload.net/images/150316/pwhoxcoe.jpg
                for (int j = 0; j < bytesSnp; j++) {
                    int[] snp_expanded = GenoCodes.gencodeMap.get(new Byte(snpBuffer[j]));
                    for(int k = 4*j, m = 0; k < 4*(j+1); k++, m++) {
                        if(k >= nIndividuals) {
                            break;
                        }
                        geno[i][k] = snp_expanded[m];
                    }
                }
            }
        } catch (FileNotFoundException e) {
            System.out.printf("File not found: %s%n", bedFilename);
        } catch (IOException e) {
            System.out.println("IO error. ");
            e.printStackTrace();
        } finally {
            bedFile.close();
        }
        if(transpose) {
            // transpose the genotype matrix
            for(int i = 0; i < nSNPtoRead; i++) {
                for(int j= 0; j < nIndividuals; j++) {
                    geno_transposed[j][i] = geno[i][j];
                }
            }
            return geno_transposed;
        } else {
            return geno;
        }
    }

    public int[][] readBed(int[] snp_vec) throws IOException {
        return readBed(snp_vec, true);
    }

    public int[][] readBed() throws IOException {
        return readBed(Utils.seq(1, nSNPs), true);
    }

    public int[][] readBed(boolean transpose) throws IOException {
        return readBed(Utils.seq(1, nSNPs), transpose);
    }

    public int[][] readBed(String[] snp_vec, boolean transpose) throws IOException {
        return readBed(getSnpIndex(snp_vec), transpose);
    }

    public int[][] readBed(String[] snp_vec) throws IOException {
        return readBed(getSnpIndex(snp_vec));
    }


    // reads the ith SNP
    public byte[] readSnpBytes(RandomAccessFile bedFile, int i) throws IOException {
        byte[] snp_buffer = new byte[bytesSnp];
        bedFile.seek((long) i * (long) bytesSnp + 3);
        bedFile.readFully(snp_buffer);
        return snp_buffer;
    }

    // reads nBytes bytes from inputStream
    public byte[] readSnpBytes(InputStream inputStream, int nBytes) throws IOException {
        byte[] buffer = new byte[nBytes];
        inputStream.read(buffer);
        return buffer;
    }

    public byte[] readSnpBytes(InputStream inputStream) throws IOException {
        return readSnpBytes(inputStream, bytesSnp);
    }


    public void copyBed(String inputFile, String outFile) throws IOException {
        FileInputStream fileInputStream = new FileInputStream(inputFile);
        BufferedInputStream bufferedInputStream = new BufferedInputStream(fileInputStream);
        FileOutputStream fileOutputStream = new FileOutputStream(outFile);
        BufferedOutputStream bufferedOutputStream = new BufferedOutputStream(fileOutputStream);
        byte[] bytesBuffer;
        try{
            // read magic bytes
            bytesBuffer = readSnpBytes(bufferedInputStream, 3);
            bufferedOutputStream.write(bytesBuffer);
            for(int i = 0; i < nSNPs; i++) {
                bytesBuffer = readSnpBytes(bufferedInputStream);
                bufferedOutputStream.write(bytesBuffer);
            }
        } finally {
            bufferedInputStream.close();
            bufferedOutputStream.close();
        }
        FileUtils.copyFile(new File(bimFilename(inputFile)), new File(bimFilename(outFile)));
        FileUtils.copyFile(new File(famFilename(inputFile)), new File(famFilename(outFile)));
    }

    public void copyBed(String outFile) throws IOException {
        copyBed(bedFilename, outFile);
    }


    public void shift(String inputFile, int nShift) throws IOException {
        shift(inputFile, nShift, GenoBytes.defaultCollpaseMatrix);
    }

    public void shift(int nShift) throws IOException {
        shift(bedFilename, nShift, GenoBytes.defaultCollpaseMatrix);
    }

    public void shift(int nShift, int[][] collapseMatrix) throws IOException {
        shift(bedFilename, nShift, collapseMatrix);
    }

    public void shift(String inputFilename, int nShift, int[][] collapseMatrix) throws IOException {
        ShiftCodes shiftCodes = new ShiftCodes(collapseMatrix);
        String outFile = Utils.shiftedFilename(inputFilename, nShift);
        String outFileBim = FilenameUtils.removeExtension(outFile) + ".bim";
        String outFileFam = FilenameUtils.removeExtension(outFile) + ".fam";
        if((new File(outFile)).exists() &&
                (new File(outFileBim).exists()) &&
                (new File(outFileFam).exists())) {
            System.out.printf("From Java, file already exists: %s", outFile);
            return;
        }

        // two input streams needed. the second with an offset
        FileInputStream fileInputStream = new FileInputStream(inputFilename);
        BufferedInputStream bufferedInputStream = new BufferedInputStream(fileInputStream);
        FileInputStream fileInputStream1 = new FileInputStream(inputFilename);
        BufferedInputStream bufferedInputStream1 = new BufferedInputStream(fileInputStream1);

        FileOutputStream fileOutputStream = new FileOutputStream(outFile);
        BufferedOutputStream bufferedOutputStream = new BufferedOutputStream(fileOutputStream);

        byte[] bytesBuffer;
        byte[] bytesBuffer1;
        byte[] bytesCollapsed;
        try{
            // read and magic bytes
            bytesBuffer = readSnpBytes(bufferedInputStream, 3);
            bufferedOutputStream.write(bytesBuffer);
            // ignore the first nShift SNPs
            readSnpBytes(bufferedInputStream1, nShift * bytesSnp + 3);
            for(int i = 0; i < nSNPs; i++) {
                bytesBuffer = readSnpBytes(bufferedInputStream);
                if(i >= nSNPs - nShift) {
                    bytesCollapsed = new byte[bytesSnp];
                    Arrays.fill(bytesCollapsed, (byte) 0b01010101);
                } else {
                    bytesBuffer1 = readSnpBytes(bufferedInputStream1);
                    bytesCollapsed = shiftCodes.lookup(bytesBuffer, bytesBuffer1);
                }
                bufferedOutputStream.write(bytesCollapsed);
            }
        } finally {
            bufferedInputStream.close();
            bufferedOutputStream.close();
        }
        String bimFilename = bimFilename(inputFilename);
        String famFilename = famFilename(inputFilename);
        String bimFilenameShifted = Utils.shiftedFilename(bimFilename, nShift);
        String famFilenameShifted = Utils.shiftedFilename(famFilename, nShift);
        try {
            Files.createSymbolicLink(Paths.get(bimFilenameShifted), Paths.get(bimFilename));
        } catch (FileAlreadyExistsException e) {
        } catch (IOException e) {
            System.err.println(e);
        }

        try{
            Files.createSymbolicLink(Paths.get(famFilenameShifted), Paths.get(famFilename));
        } catch (FileAlreadyExistsException e) {
        } catch (IOException e) {
            System.err.println(e);
        }
    }


    public String shiftedFilename(int nShift) {
        return Utils.shiftedFilename(bedFilename, nShift);
    }

    public OutputStream writeBytes(OutputStream fileStream, byte[] bytes) {
        try{
            fileStream.write(bytes);
        } catch (IOException e) {
            e.printStackTrace();
        }
        return fileStream;
    }

    public OutputStream writeBytes(OutputStream fileStream) {
        return writeBytes(fileStream, magicBytes);
    }

    public OutputStream writeBytes(String filename, byte[] bytes) throws FileNotFoundException {
        FileOutputStream fileOutputStream = new FileOutputStream(filename);
        BufferedOutputStream fileBufferedOutput = new BufferedOutputStream(fileOutputStream);
        return writeBytes(fileBufferedOutput, bytes);
    }


    public OutputStream writeBytes(String filename) throws FileNotFoundException {
        return writeBytes(filename, magicBytes);
    }

    public String getBedFilename() {
        return bedFilename;
    }

    public String getBimFilename() {
        return bimFilename;
    }

    public String getFamFilename() {
        return famFilename;
    }

    public int getBytesSnp() {
        return bytesSnp;
    }

    public int getnIndividuals() {
        return nIndividuals;
    }

    public int getnIndividualsApparent() {
        return nIndividualsApparent;
    }

    public long getBedSize() {
        return bedSize;
    }

    public int getnSNPs() {
        return nSNPs;
    }

    public byte[] getMagicBytes() {
        return magicBytes;
    }

    public static void main(String[] args) throws IOException {
        String filename = Bed.class.getResource("/test.bed").getFile();
        System.out.println(filename);

        Bed bed = new Bed(filename);

        // read in data, SNP in columns
        int[][] res;
        res = bed.readBed();
        Utils.printMat(res);
        System.out.println(bed.getSnpIndex("snp2"));
        int[] bedIndices245 = bed.getSnpIndex(new String[] {"snp2", "snp4", "snp5"});
        for(int i : bedIndices245) {
            System.out.println(i);
        }

//        String filename2 = Utils.shiftedFilename(filename, 1);
//        System.out.println(filename2);
//
//        // read in the data, SNPs in rows:
//        res = bed.readBed(false);
//        Utils.printMat(res);
//
//        bed.shift(1);
//        Bed bed2 = new Bed(filename2);
//        Utils.printMat(bed2.readBed());
//
//        String filename3 = "/Users/kaiyin/Desktop/mmp/mmp13.bed";
//        Bed bed3 = new Bed(filename3);
//        Utils.printMat(bed3.readBed(Utils.seq(1, 10)), 0, 20, 0, 10);
//        bed3.shift(1);
//        String filename3Shifted = Utils.shiftedFilename(filename3, 1);
//        Bed bed3Shifted = new Bed(filename3Shifted);
//        Utils.printMat(bed3Shifted.readBed(Utils.seq(1, 10)), 0, 20, 0, 10);
    }
}
