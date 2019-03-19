#define X 1
#define Y 2 

double minimum_coordinate(int coordinate, dmatrix_t P[], int n) {

  int i ;
  double min ;

  min = P[0].m[coordinate][1] ;

  for (i = 1 ; i < n ; i++) {
    if (P[i].m[coordinate][1] < min) {
      min = P[i].m[coordinate][1] ;
    }
  }
  return min ;
}


double maximum_coordinate(int coordinate, dmatrix_t P[], int n) {

  int i ;
  double max ;

  max = P[0].m[coordinate][1] ;

  for (i = 1 ; i < n ; i++) {
    if (P[i].m[coordinate][1] > max) {
      max = P[i].m[coordinate][1] ;
    }
  }
  return max ;
}


int maximum_intersection(int intersections[], int n) {

  int i, max ;
  
  max = intersections[0] ;
  for (i = 1 ; i < n ; i++) {
    if (intersections[i] > max) {
      max = intersections[i] ;
    }
  }
  return max ;
}


int minimum_intersection(int intersections[], int n) {

  int i, min ;
  
  min = intersections[0] ;
  for (i = 1 ; i < n ; i++) {
    if (intersections[i] < min) {
      min = intersections[i] ;
    }
  }
  return min ;
}

void XFillConvexPolygon(struct coordinate polygonLight, dmatrix_t P[], SDL_Renderer * renderer, int n) {
  
  int i, j ;
  int y, y_min, y_max, min_int, max_int ;
  double m, b ;
  int *active, *horizontal, *intersections ;

  horizontal = (int *)malloc(n*sizeof(int)) ; // Allocate horizontal segment table
  active = (int *)malloc(n*sizeof(int)) ; // Allocate active segment table
  intersections = (int *)malloc(n*sizeof(int)) ; // Allocate intersection table

  y_min = (int)minimum_coordinate(Y,P,n) ; // Determine number of scan lines
  y_max = (int)maximum_coordinate(Y,P,n) ;

  for (i = 0 ; i < n ; i++) {
   horizontal[i] = (int)P[i].m[Y][1] == (int)P[(i+1)%n].m[Y][1] ; // Find horizontal segments
  }
 
  for (y = y_min ; y <= y_max ; y++) { // For each scan line y
    for (i = 0 ; i < n ; i++) {  // Update segment table
      if (!horizontal[i]) {
        active[i] = (y >= (int)P[i].m[Y][1] && y <= (int)P[(i+1)%n].m[Y][1]) || (y <= (int)P[i].m[Y][1] && y >= (int)P[(i+1)%n].m[Y][1]) ;
      }
    }
    j = 0 ; 
    for (i = 0 ; i < n ; i++) { // find intersection x-value. The y-value is given by the scan line
      if (active[i] && !horizontal[i]) {
        if ((int)P[i].m[X][1] == (int)P[(i+1)%n].m[X][1]) { // Vertical segment
          intersections[j++] = (int)P[i].m[X][1] ; 
        }
        else {
          m = (double)((int)P[(i+1)%n].m[Y][1] - (int)P[i].m[Y][1])/(double)((int)P[(i+1)%n].m[X][1] - (int)P[i].m[X][1]) ; // Compute slope and intercept
          b = (double)((int)P[i].m[Y][1]) - m*(double)((int)P[i].m[X][1]) ;
          intersections[j++] = (int)(((double)y - b)/m) ; // Compute intersection
        }  
      }
    }
    min_int = minimum_intersection(intersections,j) ;
    max_int = maximum_intersection(intersections,j) + 1 ;
    for ( i = min_int ; i < max_int ; i++) { // Tracing from minimum to maximum intersection
      SDL_SetRenderDrawColor(renderer, (255 * polygonLight.x), (255 * polygonLight.y), (255 * polygonLight.z), 255);
      SDL_RenderDrawPoint(renderer, i, y);
    }
  }
  free(horizontal) ;
  free(active) ;
  free(intersections) ;
}
