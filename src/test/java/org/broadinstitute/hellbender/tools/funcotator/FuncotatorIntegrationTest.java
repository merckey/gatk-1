package org.broadinstitute.hellbender.tools.funcotator;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.xsv.SimpleKeyXsvFuncotationFactory;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
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
    private static final boolean doDebugTests = true;

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
        return getOutputFile( outputFormatType, "funcotator_tmp_out", outputFormatType.toString().toLowerCase() );
    }

    private static File getOutputFile(final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType,
                                      final String outfileBaseName,
                                      final String outFileExtension) {
        final File outputFile;
        if ( !doDebugTests ) {
            outputFile = createTempFile(tmpOutDir + File.separator + outfileBaseName, "." + outFileExtension);
        }
        else {
            outputFile = new File(tmpOutDir, outfileBaseName + "." + outputFormatType.toString().toLowerCase());
        }
        return outputFile;
    }

    //==================================================================================================================

    @DataProvider
    Object[][] provideForIntegrationTest() {
        return new Object[][] {
                {FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER, FuncotatorArgumentDefinitions.ReferenceVersionType.hg19, FuncotatorTestConstants.HG19_CHR3_REFERENCE_FILE_NAME, FuncotatorTestConstants.VARIANT_FILE_HG19_CHR3, FuncotatorTestConstants.PIK3CA_TRANSCRIPT, SimpleKeyXsvFuncotationFactory.XsvDataKeyType.GENE_NAME, FuncotatorArgumentDefinitions.OutputFormatType.VCF},
                {FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER, FuncotatorArgumentDefinitions.ReferenceVersionType.hg19, FuncotatorTestConstants.HG19_CHR3_REFERENCE_FILE_NAME, FuncotatorTestConstants.VARIANT_FILE_HG19_CHR3, FuncotatorTestConstants.PIK3CA_TRANSCRIPT, SimpleKeyXsvFuncotationFactory.XsvDataKeyType.TRANSCRIPT_ID, FuncotatorArgumentDefinitions.OutputFormatType.VCF},
                {FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER, FuncotatorArgumentDefinitions.ReferenceVersionType.hg19, FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME, FuncotatorTestConstants.VARIANT_FILE_HG19_CHR19, FuncotatorTestConstants.MUC16_TRANSCRIPT, SimpleKeyXsvFuncotationFactory.XsvDataKeyType.GENE_NAME, FuncotatorArgumentDefinitions.OutputFormatType.VCF},
                {FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER, FuncotatorArgumentDefinitions.ReferenceVersionType.hg19, FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME, FuncotatorTestConstants.VARIANT_FILE_HG19_CHR19, FuncotatorTestConstants.MUC16_TRANSCRIPT, SimpleKeyXsvFuncotationFactory.XsvDataKeyType.TRANSCRIPT_ID, FuncotatorArgumentDefinitions.OutputFormatType.VCF},
                {FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER, FuncotatorArgumentDefinitions.ReferenceVersionType.hg19, FuncotatorTestConstants.HG19_CHR3_REFERENCE_FILE_NAME, FuncotatorTestConstants.VARIANT_FILE_HG19_CHR3, FuncotatorTestConstants.PIK3CA_TRANSCRIPT, SimpleKeyXsvFuncotationFactory.XsvDataKeyType.GENE_NAME, FuncotatorArgumentDefinitions.OutputFormatType.MAF},
        };
    }

    //==================================================================================================================

    // This test is to make sure we don't create a bunch of temp files anywhere.
    // It will force anyone who changes the outputToTmpDir flag to make it true when they check in this test file.
    @Test
    public void metaTestEnsureTempDirs() {
        Assert.assertEquals(doDebugTests, false);
    }

    @Test(dataProvider = "provideForIntegrationTest")
    public void basicMarbleRoll(final String dataSourcesPath,
                                final FuncotatorArgumentDefinitions.ReferenceVersionType refVer,
                                final String referenceFileName,
                                final String variantFileName,
                                final String transcriptName,
                                final SimpleKeyXsvFuncotationFactory.XsvDataKeyType xsvMatchType,
                                final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType) throws IOException {

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

        for ( final FuncotatorArgumentDefinitions.OutputFormatType outFormat : Arrays.asList(FuncotatorArgumentDefinitions.OutputFormatType.VCF, FuncotatorArgumentDefinitions.OutputFormatType.MAF)) {

            final File outputFile = getOutputFile(outFormat);

            final List<String> arguments = new ArrayList<>();

            arguments.add("-" + StandardArgumentDefinitions.VARIANT_SHORT_NAME);

//        arguments.add("/Users/jonn/Development/oncotator_testing/BENCHMARK_INPUT.funcotator.vcf");
            arguments.add("/Users/jonn/Development/M2_01115161-TA1-filtered.vcf");

            arguments.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
            arguments.add("/Users/jonn/Development/references/Homo_sapiens_assembly19.fasta");

            arguments.add("--" + FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME);
            arguments.add(FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER);

            arguments.add("--" + FuncotatorArgumentDefinitions.ALLOW_HG19_GENCODE_B37_CONTIG_MATCHING_LONG_NAME);

            arguments.add("--" + FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME);
            arguments.add(FuncotatorArgumentDefinitions.ReferenceVersionType.hg19.toString());
            arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
            arguments.add(outputFile.getAbsolutePath());
            arguments.add("-" + FuncotatorArgumentDefinitions.OUTPUT_FORMAT_LONG_NAME);
            arguments.add(outFormat.toString());

            runCommandLine(arguments);
        }
    }

    @Test(dataProvider = "provideForIntegrationTest")
    public void exhaustiveArgumentTest(final String dataSourcesPath,
                                       final FuncotatorArgumentDefinitions.ReferenceVersionType refVer,
                                       final String referenceFileName,
                                       final String variantFileName,
                                       final String transcriptName,
                                       final SimpleKeyXsvFuncotationFactory.XsvDataKeyType xsvMatchType,
                                       final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType) throws IOException {

        final String outFileName = "funcotator_tmp_out_" + xsvMatchType.toString() + "_" + transcriptName + "." + outputFormatType.toString().toLowerCase();

        final File outputFile = getOutputFile( outputFormatType, outFileName.substring(0,outFileName.length()-4), outFileName.substring(outFileName.length()-3) );

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
    }
}
