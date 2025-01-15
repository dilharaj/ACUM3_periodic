void redistribute(float* t,float* a,float* b,float time,int n);
void interp_1d(float* radin,float* datin,float* radout,float* datout,int nrin,int nrout);
void phase_shift(float* u, int N, float angle);
void differentiate(float* u, float* ud, float dt, int n);
void differentiate_biased(float* u, float* ud, float dt, int n);
void make_periodic(float* fun, int n, int is,float* f_cubic);
void medfilt1(float* func, int n);
void pre_gaussfilt1(float* gauss, int n_win);
void gaussfilt1(float* func, float* gauss, int n, int n_win);
