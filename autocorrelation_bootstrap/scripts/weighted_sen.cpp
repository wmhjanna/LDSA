#include <cstdint>
#include <cmath>

static double median_counted(const double* v, const int32_t* order, const int32_t* c, int n, int64_t total) {
    const int64_t lower = (total - 1) / 2, upper = total / 2;
    int64_t cumulative = 0;
    double a = NAN;
    for (int k = 0; k < n; ++k) {
        int j = order[k];
        cumulative += c[j];
        if (cumulative > lower && std::isnan(a)) a = v[j];
        if (cumulative > upper) return (a + v[j]) / 2.0;
    }
    return NAN;
}

extern "C" void weighted_sen(int n, int m, int B,
        const int32_t* pi, const int32_t* pj, const double* slopes,
        const double* x, const double* y, const int32_t* ox, const int32_t* oy,
        const int32_t* counts, double* output) {
    for (int b = 0; b < B; ++b) {
        const int32_t* c = counts + b*n;
        int64_t total = 0, squares = 0;
        for(int j=0;j<n;++j) { total += c[j]; squares += int64_t(c[j])*c[j]; }
        const int64_t pairs = (total*total - squares)/2;
        double slope = NAN;
        if(pairs > 0) {
            const int64_t lower=(pairs-1)/2, upper=pairs/2;
            int64_t cumulative=0;
            double a=NAN;
            for(int k=0;k<m;++k) {
                cumulative += int64_t(c[pi[k]])*c[pj[k]];
                if(cumulative>lower && std::isnan(a)) a=slopes[k];
                if(cumulative>upper) { slope=(a+slopes[k])/2.0; break; }
            }
        }
        const double start = median_counted(y,oy,c,n,total)-slope*median_counted(x,ox,c,n,total);
        output[3*b]=slope;
        output[3*b+1]=start;
        output[3*b+2]=36500.0*slope/start;
    }
}
