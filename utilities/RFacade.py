import rpy2
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.vectors import FloatVector, StrVector
# r-libraries
from rpy2.robjects.packages import importr

pandas2ri.activate() # to convert pandas to R dataframe

class RFacade:
    @staticmethod
    def get_bh_pvalues(pvals):
        """
        Return a list with BH-corrected p-values
        """

        stats = importr('stats')
        p_adjust = stats.p_adjust(FloatVector(pvals), method = 'BH')
        return p_adjust

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

