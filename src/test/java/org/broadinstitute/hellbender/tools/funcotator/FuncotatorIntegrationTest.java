package org.broadinstitute.hellbender.tools.funcotator;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.xsv.SimpleKeyXsvFuncotationFactory;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * An integration test for the {@link Funcotator} tool.
 * Created by jonn on 8/29/17.
 */
public class FuncotatorIntegrationTest extends CommandLineProgramTest {

    //TODO: SUPER IMPORTANT! Add checks on the output file to make sure it's equal to an expected output file!

    // Temp directory in which to place output files.
    private static final File tmpOutDir;

    // Whether to do debug output (i.e. leave output around).
    // This should always be true when checked in.
    private static final boolean doDebugTests = false;

    static {
        if ( !doDebugTests ) {
            tmpOutDir = createTempDir("funcotatorTmpFolder");
        }
        else {
            tmpOutDir = new File("funcotatorTmpFolder" + File.separator);
            if ( !tmpOutDir.mkdirs() && !tmpOutDir.exists() ) {
                throw new GATKException("Error making output folder for test: " + FuncotatorIntegrationTest.class.getName());
            }
        }
    }

    //==================================================================================================================
    // Helper methods to create output files and maybe leave them around to debug the test output.

    private static File getOutputFile(final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType) {
        return getOutputFile( "funcotator_tmp_out", outputFormatType.toString().toLowerCase() );
    }

    private static File getOutputFile(final String outfileBaseName,
                                      final String outFileExtension) {
        final File outputFile;
        if ( !doDebugTests ) {
            outputFile = createTempFile(tmpOutDir + File.separator + outfileBaseName, "." + outFileExtension);
        }
        else {
            outputFile = new File(tmpOutDir, outfileBaseName + "." + outFileExtension);
        }
        return outputFile;
    }

    //==================================================================================================================

    @DataProvider
    Iterator<Object[]> provideForIntegrationTest() {

        final ArrayList<Object[]> testCases = new ArrayList<>();

        // TODO: Fix this set of tests!  THEY DON'T ALL PASS!  (Issue: https://github.com/broadinstitute/gatk/issues/4344)
        // VCF SNP / MNP / INDEL tests
//        for ( final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType : FuncotatorArgumentDefinitions.OutputFormatType.values() ) {
//            for ( final SimpleKeyXsvFuncotationFactory.XsvDataKeyType keyType : SimpleKeyXsvFuncotationFactory.XsvDataKeyType.values() ) {
//                testCases.add(
//                        new Object[] {
//                                FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
//                                FuncotatorArgumentDefinitions.ReferenceVersionType.hg19,
//                                FuncotatorTestConstants.HG19_CHR3_REFERENCE_FILE_NAME,
//                                FuncotatorTestConstants.PIK3CA_SNP_FILE_BASE_NAME + ".vcf",
//                                FuncotatorTestConstants.PIK3CA_TRANSCRIPT,
//                                keyType,
//                                outputFormatType,
//                                FuncotatorTestConstants.PIK3CA_SNP_FILE_BASE_NAME + ".oncotatorAnnotated." + outputFormatType.toString().toLowerCase()
//                        }
//                );
//                testCases.add(
//                        new Object[] {
//                                FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
//                                FuncotatorArgumentDefinitions.ReferenceVersionType.hg19,
//                                FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME,
//                                FuncotatorTestConstants.MUC16_MNP_FILE_BASE_NAME + ".vcf",
//                                FuncotatorTestConstants.MUC16_TRANSCRIPT,
//                                keyType,
//                                outputFormatType,
//                                FuncotatorTestConstants.MUC16_MNP_FILE_BASE_NAME + ".oncotatorAnnotated." + outputFormatType.toString().toLowerCase()
//                        }
//                );
//                testCases.add(
//                        new Object[] {
//                                FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
//                                FuncotatorArgumentDefinitions.ReferenceVersionType.hg19,
//                                FuncotatorTestConstants.HG19_CHR3_REFERENCE_FILE_NAME,
//                                FuncotatorTestConstants.PIK3CA_INDEL_FILE_BASE_NAME + ".vcf",
//                                FuncotatorTestConstants.PIK3CA_TRANSCRIPT,
//                                keyType,
//                                outputFormatType,
//                                FuncotatorTestConstants.PIK3CA_INDEL_FILE_BASE_NAME + ".oncotatorAnnotated." + outputFormatType.toString().toLowerCase()
//                        }
//                );
//            }
//        }

        // Basic Test Cases:
        for ( final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType : FuncotatorArgumentDefinitions.OutputFormatType.values() ) {
            for ( final SimpleKeyXsvFuncotationFactory.XsvDataKeyType keyType : SimpleKeyXsvFuncotationFactory.XsvDataKeyType.values() ) {
                testCases.add(
                        new Object[]{
                                FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
                                FuncotatorArgumentDefinitions.ReferenceVersionType.hg19,
                                FuncotatorTestConstants.HG19_CHR3_REFERENCE_FILE_NAME,
                                FuncotatorTestConstants.VARIANT_FILE_HG19_CHR3,
                                FuncotatorTestConstants.PIK3CA_TRANSCRIPT,
                                keyType,
                                outputFormatType,
                                null
                        }
                );
                testCases.add(
                        new Object[]{
                                FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
                                FuncotatorArgumentDefinitions.ReferenceVersionType.hg19,
                                FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME,
                                FuncotatorTestConstants.VARIANT_FILE_HG19_CHR19,
                                FuncotatorTestConstants.MUC16_TRANSCRIPT,
                                keyType,
                                outputFormatType,
                                null
                        }
                );
            }
        }

        return testCases.iterator();
    }

    //==================================================================================================================

    // This test is to make sure we don't create a bunch of temp files anywhere.
    // It will force anyone who changes the outputToTmpDir flag to make it true when they check in this test file.
    @Test
    public void metaTestEnsureTempDirs() {
        Assert.assertEquals(doDebugTests, false);
    }

    @Test(dataProvider = "provideForIntegrationTest")
    public void testFuncotatorWithoutValidatingResults(final String dataSourcesPath,
                                final FuncotatorArgumentDefinitions.ReferenceVersionType refVer,
                                final String referenceFileName,
                                final String variantFileName,
                                final String transcriptName,
                                final SimpleKeyXsvFuncotationFactory.XsvDataKeyType xsvMatchType,
                                final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType, 
                                final String expectedOutputFilePath) throws IOException {

        final File outputFile = getOutputFile(outputFormatType);

        final List<String> arguments = new ArrayList<>();

        arguments.add("--" + FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME);
        arguments.add(dataSourcesPath);
        arguments.add("--" + FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME);
        arguments.add(refVer.toString());
        arguments.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
        arguments.add(referenceFileName);
        arguments.add("-" + StandardArgumentDefinitions.VARIANT_SHORT_NAME);
        arguments.add(variantFileName);
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        arguments.add("-" + FuncotatorArgumentDefinitions.OUTPUT_FORMAT_LONG_NAME);
        arguments.add(outputFormatType.toString());

        runCommandLine(arguments);
    }

    @Test(enabled = doDebugTests)
    public void spotCheck() throws IOException {

        long startTime = 0, endTime = 0;
        final long overallStartTime = System.nanoTime();

        final List<FuncotatorArgumentDefinitions.OutputFormatType> outFormatList = new ArrayList<>();
//        outFormatList.add(FuncotatorArgumentDefinitions.OutputFormatType.VCF);
        outFormatList.add(FuncotatorArgumentDefinitions.OutputFormatType.MAF);

        for ( final FuncotatorArgumentDefinitions.OutputFormatType outFormat : outFormatList) {

            startTime = System.nanoTime();

            final File outputFile;
            if ( outFormat == FuncotatorArgumentDefinitions.OutputFormatType.VCF ) {
                outputFile = getOutputFile("funcotator_tmp_out_spot_check", outFormat.toString().toLowerCase());
            }
            else {
                outputFile = getOutputFile("funcotator_tmp_out_spot_check.maf", "tsv");
            }

            final List<String> arguments = new ArrayList<>();

            arguments.add("-" + StandardArgumentDefinitions.VARIANT_SHORT_NAME);

//        arguments.add("/Users/jonn/Development/oncotator_testing/BENCHMARK_INPUT.funcotator.vcf");
            arguments.add("/Users/jonn/Development/M2_01115161-TA1-filtered.vcf");
//            arguments.add("/Users/jonn/Development/M2_01115161-TA1-filtered.oneRecord.vcf");
//            arguments.add("/Users/jonn/Development/M2_01115161-TA1-filtered.badRecord.vcf");
//
            arguments.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
            arguments.add("/Users/jonn/Development/references/Homo_sapiens_assembly19.fasta");

            arguments.add("--" + FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME);
//            arguments.add(FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER);
            arguments.add("/Users/jonn/Development/funcotator_dataSources.v1.0.20180105");

            arguments.add("--" + FuncotatorArgumentDefinitions.ALLOW_HG19_GENCODE_B37_CONTIG_MATCHING_LONG_NAME);

            arguments.add("--" + FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME);
            arguments.add("Center:The Broad Institute");
            arguments.add("--" + FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME);
			arguments.add("dbSNP_RS:0");
            arguments.add("--" + FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME );
			arguments.add("dbSNP_Val_Status:No Value");
            arguments.add("--" + FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME );
			arguments.add("tumor_barcode:01100110011");
            arguments.add("--" + FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME );
			arguments.add("Matched_Norm_Sample_Barcode:1001");
            arguments.add("--" + FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME );
			arguments.add("Match_Norm_Seq_Allele1:No Match");
            arguments.add("--" + FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME );
			arguments.add("Match_Norm_Seq_Allele2:Abnormal");
            arguments.add("--" + FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME );
			arguments.add("Tumor_Validation_Allele1:Cancer Maybe");
            arguments.add("--" + FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME );
			arguments.add("Tumor_Validation_Allele2:Maybe Cancer");
            arguments.add("--" + FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME );
			arguments.add("Match_Norm_Validation_Allele1:Everythings fine here");
            arguments.add("--" + FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME );
			arguments.add("Match_Norm_Validation_Allele2:How are you?");
            arguments.add("--" + FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME );
			arguments.add("Verification_Status:Unverified");
            arguments.add("--" + FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME );
			arguments.add("Validation_Status:Totally Valid.  100%  Definitely.");
            arguments.add("--" + FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME );
			arguments.add("Mutation_Status:Mutated");
            arguments.add("--" + FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME );
			arguments.add("Sequencing_Phase:90 Degrees");
            arguments.add("--" + FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME );
			arguments.add("Sequence_Source:Mammal");
            arguments.add("--" + FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME );
			arguments.add("Validation_Method:Eyeballed it.");

            arguments.add("--" + FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME);
            arguments.add(FuncotatorArgumentDefinitions.ReferenceVersionType.hg19.toString());
            arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
            arguments.add(outputFile.getAbsolutePath());
            arguments.add("-" + FuncotatorArgumentDefinitions.OUTPUT_FORMAT_LONG_NAME);
            arguments.add(outFormat.toString());

            runCommandLine(arguments);

            endTime = System.nanoTime();

            System.out.println("  Elapsed Time (" + outFormat.toString() + "): " + (endTime - startTime)/1e9 + "s");
        }

        System.out.println("Total Elapsed Time: " + (endTime - overallStartTime)/1e9 + "s");
    }

    @Test(dataProvider = "provideForIntegrationTest")
    public void exhaustiveArgumentTest(final String dataSourcesPath,
                                       final FuncotatorArgumentDefinitions.ReferenceVersionType refVer,
                                       final String referenceFileName,
                                       final String variantFileName,
                                       final String transcriptName,
                                       final SimpleKeyXsvFuncotationFactory.XsvDataKeyType xsvMatchType,
                                       final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType,
                                       final String expectedOutputFilePath) throws IOException {

        final String outFileName = "funcotator_tmp_out_" + xsvMatchType.toString() + "_" + transcriptName + "." + outputFormatType.toString().toLowerCase();

        final File outputFile = getOutputFile( outFileName.substring(0,outFileName.length()-4), outFileName.substring(outFileName.length()-3) );

        final List<String> arguments = new ArrayList<>();

        // Required Args:
        arguments.add("--" + FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME);
        arguments.add(dataSourcesPath);
        arguments.add("--" + FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME);
        arguments.add(refVer.toString());
        arguments.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
        arguments.add(referenceFileName);
        arguments.add("-" + StandardArgumentDefinitions.VARIANT_SHORT_NAME);
        arguments.add(variantFileName);
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        arguments.add("-" + FuncotatorArgumentDefinitions.OUTPUT_FORMAT_LONG_NAME);
        arguments.add(outputFormatType.toString());

        // Transcript selection:
        arguments.add("--" + FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_LONG_NAME);
        arguments.add(FuncotatorArgumentDefinitions.TranscriptSelectionMode.BEST_EFFECT.toString());

        arguments.add("--" + FuncotatorArgumentDefinitions.TRANSCRIPT_LIST_LONG_NAME);
        arguments.add(transcriptName);

        // Annotation Defaults and Overrides:
        arguments.add("--" + FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME);
        arguments.add("GARBAGEDAY:SPUMONI");

        arguments.add("--" + FuncotatorArgumentDefinitions.ANNOTATION_OVERRIDES_LONG_NAME);
        arguments.add("Gencode_hugoSymbol:Freddie Mercury");

        // Run the beast:
        runCommandLine(arguments);

        // Only test for content-correctness if the output file was specified:
        if ( expectedOutputFilePath != null ) {
            // Get the expected output file:
            final File expectedOutputFile = new File(expectedOutputFilePath);

            // Make sure that the actual and expected output files are the same:
            IntegrationTestSpec.assertEqualTextFiles(outputFile, expectedOutputFile, "#");
        }
    }
}
