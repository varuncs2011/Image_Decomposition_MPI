void pgmsize(char *filename, int *nx, int *ny);
void pgmread(char *filename, void *vx, int nx, int ny);
void pgmreadchunk(char *filename, void *vx, int nx, int ny, int start_x, int start_y, int offset);
void pgmwrite(char *filename, void *vx, int nx, int ny);
