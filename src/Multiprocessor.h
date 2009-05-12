
typedef double*  vector;
typedef double** matrix;

void ourerror(char *s);
matrix newmatrix(int rows, int cols);
vector newvector(int dim);
int threader(int num,int batchsize,int thread_id,char *command,int verbose);
int setup_mqm_multi(int ntraits);
int setup_mqm_permutation(int methode);
int create_dir(int runnumber);
int delete_dir(int runnumber);
int copy_files(int runnumber);
int main(int argc, char *argv[]);
int copy_trait(int num, int num_traits);