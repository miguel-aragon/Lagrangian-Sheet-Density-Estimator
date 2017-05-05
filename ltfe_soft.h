#define MIN(a, b, c) ((((a) < (b) ? (a) : (b))) < (c) ? (((a) < (b) ? (a) : (b))) : (c))
#define MAX(a, b, c) ((((a) > (b) ? (a) : (b))) > (c) ? (((a) > (b) ? (a) : (b))) : (c))
#define LOW(a,b) ((a) <  (b) ? b : a)
#define TOP(a,b) ((a) >= (b) ? b : a)

#ifdef __cplusplus
extern "C" {
#endif

  void process_tetrahedron(float **tetra, float den, float **cube, int nx, int ny, int nz);

  void sortFour(float a, float b, float c, float d, int *inds);
  
  void render_tetrahedron_in_cube(float **tetra, float den, float **cube, int nx, int ny, int nz);

  bool isocut_line(const float isoval, const float val1, const float val2, const float *vert1, const float *vert2, float *vert_inter);
  
  int  isocut_tetrahedron(float isoval, float **tetra, float iso_ver[6][3]);
  
  void render_tetrahedron(float **tetra, float *imagen, int nx, int ny, int nz);
  
  
  void write_cube(char *filename, float **cube, int nx, int ny, int nz);
  
  void raster_triangle(const float *v1, 
		       const float *v2, 
		       const float *v3, 
		       const float den, float *buffer, const int xg, const int yg, const int q);
  
  //--- TetraDen
  void load_snapshot(char *fname, int files);
  void get_densities(int ng, int nx, int ny, int nz, float **cube);
  void save_cube(char *filename, float **cube, int ng, int nx, int ny, int nz);
  void fix_periodic_tetra(double x[4], double y[4], double z[4], double box);
  double tetra_volu_ori(double xt[4], double yt[4], double zt[4]);
  void allocate_memory_particles(void);
  void allocate_memory_cube(int _nx, int _ny, int _nz);

  void raster_triangle_interpol(const float *v1,
				const float *v2,
				const float *v3,
				const float den, float *buffer, const int xg, const int yg, const int q);



#ifdef __cplusplus
}
#endif
