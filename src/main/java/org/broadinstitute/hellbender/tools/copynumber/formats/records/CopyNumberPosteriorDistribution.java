package org.broadinstitute.hellbender.tools.copynumber.formats.records;

import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.IntegerCopyNumberState;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Map;

/**
 * A record containing a copy number posterior distribution for a single interval
 */
public class CopyNumberPosteriorDistribution {

    /**
     * Posterior distribution stores in log scale
     */
    private final Map<IntegerCopyNumberState, Double> copyNumberPosteriorDistribution;

    /**
     * Tolerated error for checking equality of doubles
     */
    private static final double ERROR_TOLERANCE = 1e-4;

    /**
     * Constructor for storing copy number posterior distribution given a mapping between copy number states
     * and their respective probabilities in log scale
     *
     * @param copyNumberPosteriorDistribution copy number posterior distribution in log scale
     */
    public CopyNumberPosteriorDistribution(final Map<IntegerCopyNumberState, Double> copyNumberPosteriorDistribution) {
        this.copyNumberPosteriorDistribution = Utils.nonNull(copyNumberPosteriorDistribution);
        final double probabilitySum = this.copyNumberPosteriorDistribution.values().stream()
                .mapToDouble(FastMath::exp).sum();
        if (MathUtils.compareDoubles(probabilitySum, 1.0, ERROR_TOLERANCE) != 0) {
            throw new IllegalArgumentException("Posterior probabilities for at at least one " +
                    "posterior record do not sum up to one");
        }
    }

    /**
     * Get the probability for a given copy number state
     */
    public double getCopyNumberPosterior(final IntegerCopyNumberState integerCopyNumberState) {
        return copyNumberPosteriorDistribution.get(integerCopyNumberState);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        final CopyNumberPosteriorDistribution that = (CopyNumberPosteriorDistribution) o;

        return copyNumberPosteriorDistribution != null
                ? copyNumberPosteriorDistribution.equals(that.copyNumberPosteriorDistribution)
                : that.copyNumberPosteriorDistribution == null;

    }

    @Override
    public int hashCode() {
        return copyNumberPosteriorDistribution != null ? copyNumberPosteriorDistribution.hashCode() : 0;
    }

    @Override
    public String toString() {
        return "CopyNumberPosteriorDistribution{" +
                "copyNumberPosteriorDistribution=" + copyNumberPosteriorDistribution +
                '}';
    }
}
