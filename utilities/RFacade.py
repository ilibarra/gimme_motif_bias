import numpy as np

class RFacade:
    @staticmethod
    def get_bh_pvalues_python(p):
        """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
        p = np.asfarray(p)
        by_descend = p.argsort()[::-1]
        by_orig = by_descend.argsort()
        steps = float(len(p)) / np.arange(len(p), 0, -1)
        q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
        return q[by_orig]

    @staticmethod
    def get_pval_asterisks(pvals, vert=False,
                           symbols=["****", "***", "**", '*'],
                           thresholds=[.0001, .001, .01, .05]):
        if vert:
            return ["\n".join(symbols[0]) if pval <= thresholds[0] else
                    ("\n".join(symbols[1]) if pval < thresholds[1] else
                     ("\n".join(symbols[2]) if pval < thresholds[2] else
                      "\n".join(symbols[3]) if pval < thresholds[3] else ""))
                    for pval in pvals]

        return [symbols[0] if pval <= thresholds[0] else
                (symbols[1]if pval < thresholds[1] else
                 (symbols[2] if pval < thresholds[2] else
                  symbols[3] if pval < thresholds[3] else ""))
                for pval in pvals]



