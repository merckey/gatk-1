package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import scala.Tuple2;

import java.util.Arrays;

/**
 * This class represents a pair of inferred genomic locations on the reference whose novel adjacency is generated
 * due to a simple SV event (in other words, a simple rearrangement between two genomic locations)
 * that is suggested by the input {@link ChimericAlignment},
 * and complications as enclosed in {@link BreakpointComplications}
 * in pinning down the locations to exact base pair resolution.
 *
 * <p>
 *     It essentially represents a bi-path "big" bubble between two reference locations.
 *     One path is the "reference path" consists of the contiguous block of bases that can be extracted from the reference,
 *     if possible (i.e. no contiguous block exists between locations from difference chromosomes).
 *     The other path is encoded with the alt haplotype sequence.
 * </p>
 */
@DefaultSerializer(NovelAdjacencyAndInferredAltHaptype.Serializer.class)
public class NovelAdjacencyAndInferredAltHaptype {

    private final SimpleInterval leftJustifiedLeftRefLoc;
    private final SimpleInterval leftJustifiedRightRefLoc;

    private final StrandSwitch strandSwitch;
    private final BreakpointComplications complication;

    private final byte[] altHaplotypeSequence;

    public NovelAdjacencyAndInferredAltHaptype(final ChimericAlignment chimericAlignment, final byte[] contigSequence,
                                               final SAMSequenceDictionary referenceDictionary) {

        strandSwitch = chimericAlignment.strandSwitch;

        try {

            final BreakpointsInference inferredClass =
                    BreakpointsInference.getInferenceClass(chimericAlignment, contigSequence, referenceDictionary);

            final Tuple2<SimpleInterval, SimpleInterval> leftJustifiedBreakpoints = inferredClass.getLeftJustifiedBreakpoints();
            leftJustifiedLeftRefLoc = leftJustifiedBreakpoints._1();
            leftJustifiedRightRefLoc = leftJustifiedBreakpoints._2();

            complication = inferredClass.getComplications();

            altHaplotypeSequence = inferredClass.getInferredAltHaplotypeSequence();

        } catch (final IllegalArgumentException iaex) { // catching IAEX specifically because it is the most likely exception thrown if there's bug, this helps quickly debugging what the problem is
            throw new GATKException("Erred when inferring breakpoint location and event type from chimeric alignment:\n" +
                    chimericAlignment.toString(), iaex);
        }
    }

    public boolean hasInsertedSequence() {
        return ! complication.getInsertedSequenceForwardStrandRep().isEmpty();
    }

    public boolean hasDuplicationAnnotation() {
        return complication.indicatesRefSeqDuplicatedOnAlt();
    }

    protected NovelAdjacencyAndInferredAltHaptype(final Kryo kryo, final Input input) {
        final String contig1 = input.readString();
        final int start1 = input.readInt();
        final int end1 = input.readInt();
        this.leftJustifiedLeftRefLoc = new SimpleInterval(contig1, start1, end1);
        final String contig2 = input.readString();
        final int start2 = input.readInt();
        final int end2 = input.readInt();
        this.leftJustifiedRightRefLoc = new SimpleInterval(contig2, start2, end2);

        this.strandSwitch = StrandSwitch.values()[input.readInt()];

        final int id = input.readInt();
        this.complication = BreakpointComplications.serialize(kryo, input, id);

        final boolean altSeqIsNull = input.readBoolean();
        if (altSeqIsNull) {
            altHaplotypeSequence = null;
        } else {
            final int arraySize = input.readInt();
            altHaplotypeSequence = new byte[arraySize];
            for (int i = 0 ; i < arraySize; ++i) {
                altHaplotypeSequence[i] = input.readByte();
            }
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        NovelAdjacencyAndInferredAltHaptype that = (NovelAdjacencyAndInferredAltHaptype) o;

        if (!leftJustifiedLeftRefLoc.equals(that.leftJustifiedLeftRefLoc)) return false;
        if (!leftJustifiedRightRefLoc.equals(that.leftJustifiedRightRefLoc)) return false;
        if (strandSwitch != that.strandSwitch) return false;
        if (!complication.equals(that.complication)) return false;
        return Arrays.equals(altHaplotypeSequence, that.altHaplotypeSequence);
    }

    @Override
    public int hashCode() {
        int result = leftJustifiedLeftRefLoc.hashCode();
        result = 31 * result + leftJustifiedRightRefLoc.hashCode();
        result = 31 * result + strandSwitch.hashCode();
        result = 31 * result + complication.hashCode();
        result = 31 * result + Arrays.hashCode(altHaplotypeSequence);
        return result;
    }

    protected void serialize(final Kryo kryo, final Output output) {
        output.writeString(leftJustifiedLeftRefLoc.getContig());
        output.writeInt(leftJustifiedLeftRefLoc.getStart());
        output.writeInt(leftJustifiedLeftRefLoc.getEnd());
        output.writeString(leftJustifiedRightRefLoc.getContig());
        output.writeInt(leftJustifiedRightRefLoc.getStart());
        output.writeInt(leftJustifiedRightRefLoc.getEnd());
        output.writeInt(strandSwitch.ordinal());

        output.writeInt(BreakpointComplications.getSerializerID(complication));
        complication.serialize(kryo, output);

        if (altHaplotypeSequence==null) {
            output.writeBoolean(true);
        } else {
            output.writeBoolean(false);
            output.writeInt(altHaplotypeSequence.length);
            for (final byte b : altHaplotypeSequence) {
                output.writeByte(b);
            }
        }
    }

    public SimpleInterval getLeftJustifiedLeftRefLoc() {
        return leftJustifiedLeftRefLoc;
    }

    public SimpleInterval getLeftJustifiedRightRefLoc() {
        return leftJustifiedRightRefLoc;
    }

    public StrandSwitch getStrandSwitch() {
        return strandSwitch;
    }

    public BreakpointComplications getComplication() {
        return complication;
    }

    public byte[] getAltHaplotypeSequence() {
        return altHaplotypeSequence;
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<NovelAdjacencyAndInferredAltHaptype> {
        @Override
        public void write(final Kryo kryo, final Output output, final NovelAdjacencyAndInferredAltHaptype novelAdjacencyAndInferredAltHaptype) {
            novelAdjacencyAndInferredAltHaptype.serialize(kryo, output);
        }

        @Override
        public NovelAdjacencyAndInferredAltHaptype read(final Kryo kryo, final Input input, final Class<NovelAdjacencyAndInferredAltHaptype> klass ) {
            return new NovelAdjacencyAndInferredAltHaptype(kryo, input);
        }
    }

    /**
     * @return Intended for use in debugging and exception message only.
     */
    @Override
    public String toString() {
        return String.format("%s\t%s\t%s\t%s", leftJustifiedLeftRefLoc.toString(), leftJustifiedRightRefLoc.toString(),
                strandSwitch.name(), complication.toString());
    }
}
