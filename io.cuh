#include <string>
inline void skipSpace(const char*& p);
inline float parseFloat(const char*& p);
void readSurfaceFileFast(const std::string& filename,
                         int& E_i,
                         int& nt,
                         float* Xham, float* Uham, float* PRham, float* Nham,
                         int XPtr_i, int NPtr_i, int UPtr_i, int PRPtr_i,
                         int nTime, int nXham, int nUham, int nNham, int nPRham,
                         int iXham, int iYham, int iZham,
                         int iUXham, int iUYham, int iUZham,
                         int iPRham, int iRHOham, int iDSham,
                         int iNXham, int iNYham, int iNZham);
void read_inputs(const char* inp_file);
void get_grid_dims(const char* filname);
void read_grid(const char* filname);
void read_q(const char* filname);
void read_cl(const char* cl_file,const char* cd_file,const char* cm_file);
void read_observers(const char* filname);
void writeTimeHistory(float* pT,float* pA,float* pP);
void read_surf(int nsurf_ham);
void read_BB_inputs(const char* inp_file_BB);
void writeBBoutput(float* BBoaspl,float* BBspl);
