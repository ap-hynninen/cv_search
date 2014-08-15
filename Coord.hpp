//
// Coordinate class
// Antti-Pekka Hynninen, antti.pekka.hynninen@nrel.gov
//

struct float3 {
  float x, y, z;
};

class Coord {

private:
  
  // Number of coordinates
  int ncoord;

  // Number of shooting points
  int nshoot;

  // Coordinates: coord[0...nshoot-1][0...ncoord-1]
  float3 **coord;

  // Residue ids
  int *resid;

  void load_coord(const char *filename);

public:

  Coord(const char *coord_filename);
  ~Coord();

};
