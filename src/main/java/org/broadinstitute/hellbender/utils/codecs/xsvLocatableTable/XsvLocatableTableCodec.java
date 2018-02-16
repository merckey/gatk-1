package org.broadinstitute.hellbender.utils.codecs.xsvLocatableTable;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.LineReader;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.readers.LineIterator;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.math.NumberUtils;
import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Codec class to read from XSV (e.g. csv, tsv, etc.) files.
 * Designed specifically with use by {@link org.broadinstitute.hellbender.tools.funcotator.Funcotator} in mind.
 *
 * Files that can be parsed by the {@link XsvLocatableTableCodec} will have a sibling configuration file of the same
 * name and the `.config` extension.  This file will contain the following keys:
 *      contig
 *      start
 *      end
 *      delimiter
 *      name
 *
 * These tables are assumed to have comment lines that start with `#` and a header that has the names for each
 * column in the table as the top row.
 *
 * Two or three columns will specify the location of each row in the data (contig, start, end; start and end can be the same
 * column).
 *
 * Created by jonn on 12/4/17.
 */
public final class XsvLocatableTableCodec extends AsciiFeatureCodec<XsvTableFeature> {

    private static final Logger logger = LogManager.getLogger(XsvLocatableTableCodec.class);

    //==================================================================================================================
    // Public Static Members:

    public static final String COMMENT_DELIMITER = "#";
    public static final String SAM_HEADER_DELIMITER = "@";

    public static final String CONFIG_FILE_CONTIG_COLUMN_KEY = "contig_column";
    public static final String CONFIG_FILE_START_COLUMN_KEY = "start_column";
    public static final String CONFIG_FILE_END_COLUMN_KEY = "end_column";
    public static final String CONFIG_FILE_DELIMITER_KEY = "xsv_delimiter";
    public static final String CONFIG_FILE_DATA_SOURCE_NAME_KEY = "name";

    //==================================================================================================================
    // Private Static Members:


    private static final String CONFIG_FILE_EXTENSION = ".config";

    //==================================================================================================================
    // Private Members:

    /** Column name (or index) from which to get the contig string for each entry.  As specified in the input.*/
    private String inputContigColumn;

    /** Column name (or index) from which to get the start position for each entry.  As specified in the input. */
    private String inputStartColumn;

    /** Column name (or index) from which to get the end position for each entry.  As specified in the input. */
    private String inputEndColumn;

    /** Column name from which to get the contig string for each entry. */
    private String finalContigColumn;

    /** Column name from which to get the start position for each entry. */
    private String finalStartColumn;

    /** Column name from which to get the end position for each entry. */
    private String finalEndColumn;

    /** Delimiter for entries in this XSV Table. */
    private String delimiter;

    /** The name of the data source that is associated with this {@link XsvLocatableTableCodec}. */
    private String dataSourceName;

    /** The XSV Table Header */
    private List<String> header;

    /** The XSV Table Header */
    private Map<String, Integer> headerToIndex;

    /** The locatable columns, once determined, in contig, start, end order. */
    private List<String> locatableColumns;

    /** The current position in the file that is being read. */
    private long currentLine = 0;

    /** Config file to use instead of a sibling config file.  Null if not using an override.*/
    private Path overrideConfigFile = null;

    /** Comments, if any */
    private List<String> comments = new ArrayList<>();

    /** SAM header as strings. */
    private List<String> samFileHeaderAsStrings = new ArrayList<>();

    //==================================================================================================================
    // Constructors:

    public XsvLocatableTableCodec() {
        super(XsvTableFeature.class);
    }

    /** Constructor for when a configuration file is specified instead of using a sibling config file.
     *
     * This cannot be used with auto decoding.
     *
     * @param overrideConfigFile {@link Path} to the file to use as a configuration file for the given file.
     */
    public XsvLocatableTableCodec(final Path overrideConfigFile) {
        super(XsvTableFeature.class);
        this.overrideConfigFile = overrideConfigFile;
    }

    //==================================================================================================================
    // Override Methods:

    @Override
    public boolean canDecode(final String path) {

        // Get the paths to our file and the config file:
        final Path inputFilePath = IOUtils.getPath(path);
        final Path configFilePath = (overrideConfigFile != null ?
                overrideConfigFile : getConfigFilePath(inputFilePath));

        // Check that our files are good for eating... I mean reading...
        if ( validateInputDataFile(inputFilePath) && validateInputDataFile(configFilePath) ) {

            // Get our metadata and set up our internals so we can read from this file:
            readMetadataFromConfigFile(configFilePath);
            return true;
        }
        else {
            return false;
        }
    }

    @Override
    public XsvTableFeature decode(final String s) {

        // Increment our line counter:
        ++currentLine;

        if (s.startsWith(COMMENT_DELIMITER)) {
            return null;
        }

        if (s.startsWith(SAM_HEADER_DELIMITER)) {
            return null;
        }

        final List<String> split = new ArrayList<>(Arrays.asList(s.split(delimiter)));
        if (split.size() < 1) {
            throw new UserException.BadInput("XSV file has a line with no delimiter at line number: " + currentLine);
        }
        else if ( split.size() < header.size() ) {
            logger.warn("WARNING: Line " + currentLine + " does not have the same number of fields as header!  Padding with empty fields to end...");
            while (split.size() < header.size() ) {
                split.add("");
            }
        }
        else if ( split.size() > header.size() ) {
            logger.warn("WARNING: Line " + currentLine + " does not have the same number of fields as header!  Truncating fields from end...");
            while (split.size() > header.size() ) {
                split.remove( split.size() - 1 );
            }
        }

        return new XsvTableFeature(headerToIndex.get(finalContigColumn), headerToIndex.get(finalStartColumn),
                headerToIndex.get(finalEndColumn), header, split, dataSourceName);
    }

    /**
     * Dev note:  We also determine the actual locatable columns here.
     *
     * @param reader
     * @return
     */
    @Override
    public List<String> readActualHeader(final LineIterator reader) {
        // All leading lines with comments / header info are headers:
        while ( reader.hasNext() ) {

            final String line = reader.next();
            ++currentLine;

            // Ignore commented out lines:
            if ( !line.startsWith(COMMENT_DELIMITER) && !line.startsWith(SAM_HEADER_DELIMITER)) {

                // The first non-commented line is the column header.
                // Add the data source name to teh start of each header row,
                // then add those rows to the header object.
                header = Arrays.stream(line.split(delimiter))
                        .map(x -> determinePrefixForHeader() + x)
                        .collect(Collectors.toCollection(ArrayList::new));
                headerToIndex = IntStream.range(0, header.size()).boxed()
                        .collect(Collectors.toMap(i-> header.get(i), Function.identity()));
                finalContigColumn = NumberUtils.isNumber(inputContigColumn) ? header.get(Integer.valueOf(inputContigColumn))
                        : determinePrefixForHeader() + inputContigColumn;
                finalStartColumn = StringUtils.isNumeric(inputStartColumn) ? header.get(Integer.valueOf(inputStartColumn))
                        : determinePrefixForHeader() + inputStartColumn;
                finalEndColumn = StringUtils.isNumeric(inputEndColumn) ? header.get(Integer.valueOf(inputEndColumn))
                        : determinePrefixForHeader() + inputEndColumn;

                locatableColumns = Lists.newArrayList(finalContigColumn, finalStartColumn, finalEndColumn);

                return header;
            }

            if (line.startsWith(COMMENT_DELIMITER)) {
                comments.add(line);
            }

            if (line.startsWith(SAM_HEADER_DELIMITER)) {
                samFileHeaderAsStrings.add(line);
            }
        }

        throw new UserException.BadInput("Given file is malformed - does not contain a header!");
    }

    private String determinePrefixForHeader() {
        return (StringUtils.isEmpty(dataSourceName) ? "" : dataSourceName + "_");
    }

    //==================================================================================================================
    // Static Methods:

    /**
     * Get the properties from the given {@code configFilePath}, validate that all required properties are present,
     * and return the property map.
     * @param configFilePath {@link Path} to the configuration file.
     * @return The {@link Properties} as contained in the given {@code configFilePath}.
     */
    public static Properties getAndValidateConfigFileContents(final Path configFilePath) {

        Utils.nonNull(configFilePath);

        // Read in the contents of the config file:
        final Properties configFileContents = new Properties();
        try ( final InputStream inputStream = Files.newInputStream(configFilePath, StandardOpenOption.READ) ) {
            configFileContents.load(inputStream);
        }
        catch (final Exception ex) {
            throw new UserException.BadInput("Unable to read from XSV config file: " + configFilePath.toUri().toString(), ex);
        }

        // Validate that it has the right keys:
        assertConfigPropertiesContainsKey(configFileContents, CONFIG_FILE_CONTIG_COLUMN_KEY, configFilePath);
        assertConfigPropertiesContainsKey(configFileContents, CONFIG_FILE_START_COLUMN_KEY, configFilePath);
        assertConfigPropertiesContainsKey(configFileContents, CONFIG_FILE_END_COLUMN_KEY, configFilePath);
        assertConfigPropertiesContainsKey(configFileContents, CONFIG_FILE_DELIMITER_KEY, configFilePath);
        assertConfigPropertiesContainsKey(configFileContents, CONFIG_FILE_DATA_SOURCE_NAME_KEY, configFilePath);

        return configFileContents;
    }

    /**
     * Gets the path to the corresponding configuration file for the given {@code inputFilePath}.
     * The resulting path may or may not exist.
     * @param inputFilePath The data file {@link Path} from which to construct the path to the configuration file.
     * @return The {@link Path} for the configuration file associated with {@code inputFilePath}.
     */
    public static Path getConfigFilePath(final Path inputFilePath) {
        // Check for a sibling config file with the same name, .config as extension
        final String configFilePath = IOUtils.replaceExtension( inputFilePath.toString(), CONFIG_FILE_EXTENSION );
        return Paths.get(configFilePath);
    }

    /**
     * Ensures that the given {@link Properties} contain the given key.
     * @param configProperties The {@link Properties} in which to look for the given key.
     * @param key The value to find in the given {@link Properties}.
     * @param configFilePath The {@link Path} for the config file from which {@link Properties} were derived.  Used for printing output only.
     */
    private static void assertConfigPropertiesContainsKey(final Properties configProperties, final String key, final Path configFilePath) {
        if ( !configProperties.stringPropertyNames().contains(key) ) {
            throw new UserException.BadInput("Config file for datasource (" + configFilePath.toUri().toString() + ") does not contain required key: " + key);
        }
    }

    //==================================================================================================================
    // Instance Methods:

    /**
     * Asserts that the given {@code filePath} is a valid file from which to read.
     * @param filePath The {@link Path} to the data file to validate.
     * @return {@code true} if the given {@code filePath} is valid; {@code false} otherwise.
     */
    private boolean validateInputDataFile(final Path filePath) {
        return Files.exists(filePath) && Files.isReadable(filePath) && !Files.isDirectory(filePath);
    }

    /**
     * Reads the metadata required for parsing from the given {@code configFilePath}.
     * @param configFilePath {@link Path} to the configuration file from which to read in and setup metadata values.
     */
    private void readMetadataFromConfigFile(final Path configFilePath) {

        final Properties configProperties = getAndValidateConfigFileContents(configFilePath);

        // Get the properties and remove the leading/trailing whitespace if there is any:
        inputContigColumn = configProperties.getProperty(CONFIG_FILE_CONTIG_COLUMN_KEY).replaceAll("^\\s+", "").replaceAll("\\s+$", "");
        inputStartColumn = configProperties.getProperty(CONFIG_FILE_START_COLUMN_KEY).replaceAll("^\\s+", "").replaceAll("\\s+$", "");
        inputEndColumn = configProperties.getProperty(CONFIG_FILE_END_COLUMN_KEY).replaceAll("^\\s+", "").replaceAll("\\s+$", "");
        dataSourceName = configProperties.getProperty(CONFIG_FILE_DATA_SOURCE_NAME_KEY).replaceAll("^\\s+", "").replaceAll("\\s+$", "");

        // Get the delimiter - we do NOT remove whitespace here on purpose:
        delimiter      = configProperties.getProperty(CONFIG_FILE_DELIMITER_KEY);

        // Process delimiter just in case it is a tab escape character:
        if ( delimiter.equals("\\t") ) {
            delimiter = "\t";
        }
    }

    /** TODO: Tests of this class and overridden config files.
     * Creates a copy of the internal comments upon each invocation.
     *
     * @return an immutable list of all of the comment lines in the xsv
     */
    public ImmutableList<String> getComments() {
        return ImmutableList.copyOf(comments);
    }

    /**
     * @return copy of the sam file header created from the input file.  {@code null} is possible
     */
    public SAMFileHeader createSamFileHeader() {
        final LineReader reader = BufferedLineReader.fromString(StringUtils.join(samFileHeaderAsStrings, "\n"));
        final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
        return codec.decode(reader, null);
    }

    /**
     * Get the header from this {@link XsvLocatableTableCodec} without the columns that contain location information.
     * Specifically the columns specified by the following fields are not included:
     *  {@link XsvLocatableTableCodec#inputContigColumn}
     *  {@link XsvLocatableTableCodec#inputStartColumn}
     *  {@link XsvLocatableTableCodec#inputEndColumn}
     * @return The header for this {@link XsvLocatableTableCodec} without location columns.
     */
    public List<String> getHeaderWithoutLocationColumns() {
        return header.stream().filter(h -> !locatableColumns.contains(h))
                .collect(Collectors.toList());
    }

    public String getFinalContigColumn() {
        return finalContigColumn;
    }

    public String getFinalStartColumn() {
        return finalStartColumn;
    }

    public String getFinalEndColumn() {
        return finalEndColumn;
    }

    /**
     * Throw an exception if the given column name cannot be used for one of the locatable columns.
     * @param columnName candidate column name for one of the locatable fields (contig, start, or end)
     */
    public static void validateLocatableColumnName(String columnName) {
        Utils.validateArg(!StringUtils.isEmpty(columnName), "column header is blank.");
        Utils.validateArg(!org.apache.commons.lang.math.NumberUtils.isNumber(columnName), "column header cannot be a number: " + columnName);
    }

    //==================================================================================================================
    // Helper Data Types:
}
