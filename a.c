/*
 * CS3388 Assignment 3
 * Author: Ben Heffernan
 * Date of Creation: March 6th, 2019
 * Purpose: The main file of the program. It includes the main logic for creating the polygon mesh
 * for spheres, the lighting model using the RGB colour scheme, functionality for rendering
 * polygons based on depth and the polygon filling algorithm. The logic to create cone and
 * torus meshes was not implemented.
 */

#include <stdio.h>
#include <math.h>
#include <SDL2/SDL.h>

#include "constants.h"
#include "utils.h"

#include "camera.c"
#include "XFillPolygon.c"

/*
 * calculateDistanceFromCamera(struct coordinate pointA, struct coordinate pointB)
 * Input:
 * - pointA = Center point of a polygon
 * - pointB = Camera coordinates
 * Output: Returns an array containing all of the vertices of the polygons of the sphere polygon mesh as coordinate structs
 * Description: Calculates and places latitude and longitude lines on the sphere to create a polygon mesh of
 * squares and triangles at the poles. It stores the vertices of these polygons in a list and returns it.
*/

float calculateDistanceFromCamera(struct coordinate pointA, struct coordinate pointB)
{
    return sqrt(pow((pointA.x - pointB.x), 2.0) + pow((pointA.y - pointB.y), 2.0) + pow((pointA.z - pointB.z), 2.0));
}

/*
 * generateVertexList(point pt, float radius, int nLatitude, int nLongitude)
 * Input:
 * - pt = Center point of the sphere
 * - radius = Radius of the sphere
 * - nLatitude = Number of latitude lines to divide up the sphere horizontally for polygon creation
 * - nLongitude = Number of longitude lines to divide up the sphere vertically for polygon creation
 * Output: Returns an array containing all of the vertices of the polygons of the sphere polygon mesh as coordinate structs
 * Description: Calculates and places latitude and longitude lines on the sphere to create a polygon mesh of
 * squares and triangles at the poles. It stores the vertices of these polygons in a list and returns it.
*/

coordinateArray generateVertexList(point pt, float radius, int nLatitude, int nLongitude)
{
    int p, s, i, j;
    float x, y, z, out;
    int nPitch = nLongitude + 1;
    int numVertices = 0;

    coordinateArray sphereVertices;
    initCoordinateArray(&sphereVertices, 10);
  
    // Set up variables needed to calculate latitude and longitude lines
    float pitchInc = (180. / (float)nPitch) * DEGS_TO_RAD;
    float rotInc   = (360. / (float)nLatitude) * DEGS_TO_RAD;

    // Create the vertices for the top and bottom pole of the sphere
    struct coordinate topVertex;
    topVertex.x = pt.x;
    topVertex.y = pt.y + radius;
    topVertex.z = pt.z;

    struct coordinate bottomVertex;
    bottomVertex.x = pt.x;
    bottomVertex.y = pt.y - radius;
    bottomVertex.z = pt.z;

    insertCoordinateArray(&sphereVertices, topVertex);
    insertCoordinateArray(&sphereVertices, bottomVertex);

    numVertices = numVertices + 2;

    // Record the first vertex index for intermediate vertices
    int fVert = numVertices;

    // Loops to generate all the remaining vertices
    for (p = 1; p < nPitch; p++)
    {
        out = radius * sin((float)p * pitchInc);

        if (out < 0)
        {
            out = -out;
        }

        y = radius * cos(p * pitchInc);

        for (s = 0; s < nLatitude; s++)
        {
            x = out * cos(s * rotInc);
            z = out * sin(s * rotInc);

            struct coordinate tempVertex;
            tempVertex.x = x + pt.x;
            tempVertex.y = y + pt.y;
            tempVertex.z = z + pt.z;

            insertCoordinateArray(&sphereVertices, tempVertex);

            numVertices++;
        }
    }

    return sphereVertices;
}

/*
 * generateFaceList(int nLatitude, int nLongitude)
 * Input:
 * - nLatitude = Number of latitude lines to divide up the sphere horizontally for polygon creation
 * - nLongitude = Number of longitude lines to divide up the sphere vertically for polygon creation
 * Output: Returns an array containing all of the faces of the polygons of the sphere polygon mesh as arrays of their vertices
 * Description: Generates the face for each polygon by and stores the polygon as a list of it's vertices and then
 * returns this list of the polygons.
*/

faceArray generateFaceList(int nLatitude, int nLongitude)
{
    int p, i, j, s;
    int nPitch = nLongitude + 1;
    int fVert = 2;

    faceArray tempFaceArray;
    initFaceArray(&tempFaceArray, 50);

    // Generates the square polygon faces
    for (p = 1; p < nPitch - 1; p++)
    {
        for (s = 0; s < nLatitude; s++)
        {
            i = p * nLatitude + s;
            j = (s == nLatitude - 1) ? i - nLatitude : i;

            vertexArray tempVertexArray;
            initVertexArray(&tempVertexArray, 4);

            insertVertexArray(&tempVertexArray, (i + 1 - nLatitude) + fVert);
            insertVertexArray(&tempVertexArray, (j + 2 - nLatitude) + fVert);
            insertVertexArray(&tempVertexArray, (j + 2) + fVert);
            insertVertexArray(&tempVertexArray, (i + 1) + fVert);

            insertFaceArray(&tempFaceArray, tempVertexArray);
        }
    }

    //Generates the triangluar polygons needed to render the top and bottom of the sphere

    int offLastVerts  = fVert + (nLatitude * (nLongitude - 1));

    for(s = 0; s < nLatitude; s++)
    {
        j = (s == nLatitude - 1) ? -1 : s;

        // Get the faces for the top of the sphere
        vertexArray tempVertexArray1;
        initVertexArray(&tempVertexArray1, 3);

        insertVertexArray(&tempVertexArray1, fVert - 1);
        insertVertexArray(&tempVertexArray1, (j + 2) + fVert);
        insertVertexArray(&tempVertexArray1, (s + 1) + fVert);

        // Get the faces for the bottom of the sphere
        vertexArray tempVertexArray2;
        initVertexArray(&tempVertexArray2, 3);

        insertVertexArray(&tempVertexArray2, fVert);
        insertVertexArray(&tempVertexArray2, (s + 1) + offLastVerts);
        insertVertexArray(&tempVertexArray2, (j + 1) + offLastVerts);

        insertFaceArray(&tempFaceArray, tempVertexArray1);
        insertFaceArray(&tempFaceArray, tempVertexArray2);
    }

    return tempFaceArray;
}

/*
 * calculateTriangleCenterPoint(struct coordinate A, struct coordinate B, struct coordinate C)
 * Input:
 * - A = Coordinate of first vertex of polygon
 * - B = Coordinate of second vertex of polygon
 * - C = Coordinate of third vertex of polygon
 * Output: Returns the coordinates of the center point of a triangular polygon
 * Description: Uses the three vertices of a triangular polygon to calculate its center point
*/

struct coordinate calculateTriangleCenterPoint(struct coordinate A, struct coordinate B, struct coordinate C)
{
    // Use points A, B, and C to calculate the coordinates of the center point
    struct coordinate centerPoint;
    centerPoint.x = (A.x + B.x + C.x) / 3;
    centerPoint.y = (A.y + B.y + C.y) / 3;
    centerPoint.z = (A.z + B.z + C.z) / 3;

    return centerPoint;
}

/*
 * calculateSquareCenterPoint(struct coordinate A, struct coordinate B, struct coordinate C, struct coordinate D)
 * Input:
 * - A = Coordinate of first vertex of polygon
 * - B = Coordinate of second vertex of polygon
 * - C = Coordinate of third vertex of polygon
 * - D = Coordinate of fourth vertex of polygon
 * Output: Returns the coordinates of the center point of a rectangular polygon
 * Description: Using a ton of if statements this function figures out which two of the four
 * vertices of the polygon are furthest from each other in the x direction, in the y direction,
 * and in the z direction. It then uses the distance between the two points for each direction
 * to calculate the overall center point of the square.
*/

struct coordinate calculateSquareCenterPoint(struct coordinate A, struct coordinate B, struct coordinate C, struct coordinate D)
{
    struct coordinate centerPoint;

    int tempMaxX = fabs(A.x - B.x);
    int tempMaxY = fabs(A.y - B.y);
    int tempMaxZ = fabs(A.z - B.z);

    int ABx_bool = 1;
    int ACx_bool = 0;
    int ADx_bool = 0;
    int BCx_bool = 0;
    int BDx_bool = 0;
    int CDx_bool = 0;
    
    int ABy_bool = 1;
    int ACy_bool = 0;
    int ADy_bool = 0;
    int BCy_bool = 0;
    int BDy_bool = 0;
    int CDy_bool = 0;

    int ABz_bool = 1;
    int ACz_bool = 0;
    int ADz_bool = 0;
    int BCz_bool = 0;
    int BDz_bool = 0;
    int CDz_bool = 0;

    // Find opposite x points
    if (fabs(A.x - C.x) > tempMaxX)
    {
        tempMaxX = fabs(A.x - C.x);
        ACx_bool = 1;
        ABx_bool = 0;
    }
    if (fabs(A.x - D.x) > tempMaxX)
    {
        tempMaxX = fabs(A.x - D.x);
        ADx_bool = 1;
        ABx_bool = 0;
    }
    if (fabs(B.x - C.x) > tempMaxX)
    {
        tempMaxX = fabs(B.x - C.x);
        BCx_bool = 1;
        ABx_bool = 0;
    }
    if (fabs(B.x - D.x) > tempMaxX)
    {
        tempMaxX = fabs(B.x - D.x);
        BDx_bool = 1;
        ABx_bool = 0;
    }
    if (fabs(C.x - D.x) > tempMaxX)
    {
        tempMaxX = fabs(C.x - D.x);
        CDx_bool = 1;
        ABx_bool = 0;
    }

    // Find opposite y points
    if (fabs(A.y - C.y) > tempMaxY)
    {
        tempMaxY = fabs(A.y - C.y);
        ACy_bool = 1;
        ABy_bool = 0;
    }
    if (fabs(A.y - D.y) > tempMaxY)
    {
        tempMaxY = fabs(A.y - D.y);
        ADy_bool = 1;
        ABy_bool = 0;
    }
    if (fabs(B.y - C.y) > tempMaxY)
    {
        tempMaxY = fabs(B.y - C.y);
        BCy_bool = 1;
        ABy_bool = 0;
    }
    if (fabs(B.y - D.y) > tempMaxY)
    {
        tempMaxY = fabs(B.y - D.y);
        BDy_bool = 1;
        ABy_bool = 0;
    }
    if (fabs(C.y - D.y) > tempMaxY)
    {
        tempMaxY = fabs(C.y - D.y);
        CDy_bool = 1;
        ABy_bool = 0;
    }

    // Find opposite z points
    if (fabs(A.z - C.z) > tempMaxZ)
    {
        tempMaxZ = fabs(A.z - C.z);
        ACz_bool = 1;
        ABz_bool = 0;
    }
    if (fabs(A.z - D.z) > tempMaxZ)
    {
        tempMaxZ = fabs(A.z - D.z);
        ADz_bool = 1;
        ABz_bool = 0;
    }
    if (fabs(B.z - C.z) > tempMaxZ)
    {
        tempMaxZ = fabs(B.z - C.z);
        BCz_bool = 1;
        ABz_bool = 0;
    }
    if (fabs(B.z - D.z) > tempMaxZ)
    {
        tempMaxZ = fabs(B.z - D.z);
        BDz_bool = 1;
        ABz_bool = 0;
    }
    if (fabs(C.z - D.z) > tempMaxZ)
    {
        tempMaxZ = fabs(C.z - D.z);
        CDz_bool = 1;
        ABz_bool = 0;
    }

    // Calculate x coordinate
    if (ABx_bool == 1)
    {
        if (A.x > B.x)
        {
            centerPoint.x = A.x - (tempMaxX / 2);
        }
        else if (A.x < B.x)
        {
            centerPoint.x = A.x + (tempMaxX / 2);
        }
        else
        {
            centerPoint.x = A.x;
        }
    }
    else if (ACx_bool == 1)
    {
        if (A.x > C.x)
        {
            centerPoint.x = A.x - (tempMaxX / 2);
        }
        else if (A.x < C.x)
        {
            centerPoint.x = A.x + (tempMaxX / 2);
        }
        else
        {
            centerPoint.x = A.x;
        }
    }
    else if (ADx_bool == 1)
    {
        if (A.x > D.x)
        {
            centerPoint.x = A.x - (tempMaxX / 2);
        }
        else if (A.x < D.x)
        {
            centerPoint.x = A.x + (tempMaxX / 2);
        }
        else
        {
            centerPoint.x = A.x;
        }
    }
    else if (BCx_bool == 1)
    {
        if (B.x > C.x)
        {
            centerPoint.x = B.x - (tempMaxX / 2);
        }
        else if (B.x < C.x)
        {
            centerPoint.x = B.x + (tempMaxX / 2);
        }
        else
        {
            centerPoint.x = B.x;
        }
    }
    else if (BDx_bool == 1)
    {
        if (B.x > D.x)
        {
            centerPoint.x = B.x - (tempMaxX / 2);
        }
        else if (B.x < D.x)
        {
            centerPoint.x = B.x + (tempMaxX / 2);
        }
        else
        {
            centerPoint.x = B.x;
        }
    }
    else if (CDx_bool == 1)
    {
        if (C.x > D.x)
        {
            centerPoint.x = C.x - (tempMaxX / 2);
        }
        else if (C.x < D.x)
        {
            centerPoint.x = C.x + (tempMaxX / 2);
        }
        else
        {
            centerPoint.x = C.x;
        }
    }

    // Calculate y coordinate
    if (ABy_bool == 1)
    {
        if (A.y > B.y)
        {
            centerPoint.y = A.y - (tempMaxY / 2);
        }
        else if (A.y < B.y)
        {
            centerPoint.y = A.y + (tempMaxY / 2);
        }
        else
        {
            centerPoint.y = A.y;
        }
    }
    else if (ACy_bool == 1)
    {
        if (A.y > C.y)
        {
            centerPoint.y = A.y - (tempMaxY / 2);
        }
        else if (A.y < C.y)
        {
            centerPoint.y = A.y + (tempMaxY / 2);
        }
        else
        {
            centerPoint.y = A.y;
        }
    }
    else if (ADy_bool == 1)
    {
        if (A.y > D.y)
        {
            centerPoint.y = A.y - (tempMaxY / 2);
        }
        else if (A.y < D.y)
        {
            centerPoint.y = A.y + (tempMaxY / 2);
        }
        else
        {
            centerPoint.y = A.y;
        }
    }
    else if (BCy_bool == 1)
    {
        if (B.y > C.y)
        {
            centerPoint.y = B.y - (tempMaxY / 2);
        }
        else if (B.y < C.y)
        {
            centerPoint.y = B.y + (tempMaxY / 2);
        }
        else
        {
            centerPoint.y = B.y;
        }
    }
    else if (BDy_bool == 1)
    {
        if (B.y > D.y)
        {
            centerPoint.y = B.y - (tempMaxY / 2);
        }
        else if (B.y < C.y)
        {
            centerPoint.y = B.y + (tempMaxY / 2);
        }
        else
        {
            centerPoint.y = B.y;
        }
    }
    else if (CDy_bool == 1)
    {
        if (C.y > D.y)
        {
            centerPoint.y = C.y - (tempMaxY / 2);
        }
        else if (C.y < D.y)
        {
            centerPoint.y = C.y + (tempMaxY / 2);
        }
        else
        {
            centerPoint.y = C.y;
        }
    }

    // Calculate z coordinate
    if (ABz_bool == 1)
    {
        if (A.z > B.z)
        {
            centerPoint.z = A.z - (tempMaxZ / 2);
        }
        else if (A.z < B.z)
        {
            centerPoint.z = A.z + (tempMaxZ / 2);
        }
        else
        {
            centerPoint.z = A.z;
        }
    }
    else if (ACz_bool == 1)
    {
        if (A.z > C.z)
        {
            centerPoint.z = A.z - (tempMaxZ / 2);
        }
        else if (A.z < C.z)
        {
            centerPoint.z = A.z + (tempMaxZ / 2);
        }
        else
        {
            centerPoint.z = A.z;
        }
    }
    else if (ADz_bool == 1)
    {
        if (A.z > D.z)
        {
            centerPoint.z = A.z - (tempMaxZ / 2);
        }
        else if (A.z < D.z)
        {
            centerPoint.z = A.z + (tempMaxZ / 2);
        }
        else
        {
            centerPoint.z = A.z;
        }
    }
    else if (BCz_bool == 1)
    {
        if (B.z > C.z)
        {
            centerPoint.z = B.z - (tempMaxZ / 2);
        }
        else if (B.z < C.z)
        {
            centerPoint.z = B.z + (tempMaxZ / 2);
        }
        else
        {
            centerPoint.z = B.z;
        }
    }
    else if (BDz_bool == 1)
    {
        if (B.z > D.z)
        {
            centerPoint.z = B.z - (tempMaxZ / 2);
        }
        else if (B.z < D.z)
        {
            centerPoint.z = B.z + (tempMaxZ / 2);
        }
        else
        {
            centerPoint.z = B.z;
        }
    }
    else if (CDz_bool == 1)
    {
        if (C.z > D.z)
        {
            centerPoint.z = C.z - (tempMaxZ / 2);
        }
        else if (C.z < D.z)
        {
            centerPoint.z = C.z + (tempMaxZ / 2);
        }
        else
        {
            centerPoint.z = C.z;
        }
    }

    return centerPoint;
}

/*
 * generatePolygonCenterPointList(coordinateArray vertices, faceArray faces)
 * Input:
 * - vertices = A list of all the vertices of the polygons making up the polygon mesh for the sphere
 * - faces = A list of the faces of the polygon mesh as lists of the index of each vertex of the polygon in the vertices list
 * Output: Returns an array containing the center point for each polygon of the sphere polygon mesh
 * Description: Uses the calculateTriangleCenterPoint and calculateSquareCenterPoint functions
 * to calculate the center point for every vertex in the polygon mesh and returns them all in an array.
*/

coordinateArray generatePolygonCenterPointList(coordinateArray vertices, faceArray faces)
{
    // Array to hold the calculated center points
    coordinateArray centerPoints;
    initCoordinateArray(&centerPoints, 10);

    // Loop to get the vertex indices from faces, then to get the vertex coordinates matching each index from vertices
    for (int i = 0; i < faces.used; i++)
    {
        struct coordinate tempCenterPoint;

        if (faces.array[i].used == 3)
        {
            tempCenterPoint = calculateTriangleCenterPoint(
                vertices.array[faces.array[i].array[0]],
                vertices.array[faces.array[i].array[1]],
                vertices.array[faces.array[i].array[2]]
            );
        }
        else if (faces.array[i].used == 4)
        {
            tempCenterPoint = calculateSquareCenterPoint(
                vertices.array[faces.array[i].array[0]],
                vertices.array[faces.array[i].array[1]],
                vertices.array[faces.array[i].array[2]],
                vertices.array[faces.array[i].array[3]]
            );
        }

        struct coordinate cameraCoordinates;
        cameraCoordinates.x = Ex;
        cameraCoordinates.y = Ey;
        cameraCoordinates.z = Ez;

        tempCenterPoint.distanceFromCamera = calculateDistanceFromCamera(tempCenterPoint, cameraCoordinates);
        faces.array[i].distanceFromCamera = tempCenterPoint.distanceFromCamera;

        insertCoordinateArray(&centerPoints, tempCenterPoint);
    }

    return centerPoints;
}

/*
 * calculateFaceNormal(struct coordinate A, struct coordinate B, struct coordinate C)
 * Input:
 * - A = The first vertex of the polygon
 * - B = The second vertex of the polygon
 * - C = The third vertex of the polygon
 * Output: Returns a coordinate representing the vector of the normal of the polygon
 * Description: Uses the three points to create two vectors travelling along the edges of
 * the vertex and takes the dot product of these vectors to calculate the normal of the
 * polygon.
*/

struct coordinate calculateFaceNormal(struct coordinate A, struct coordinate B, struct coordinate C)
{
    // Use points A and B to get vector 1
    struct coordinate AB;
    AB.x = B.x - A.x;
    AB.y = B.y - A.y;
    AB.z = B.z - A.z;

    // Use points B and C to get vector 2
    struct coordinate BC;
    BC.x = C.x - B.x;
    BC.y = C.y - B.y;
    BC.z = C.z - B.z;

    // Calculate dot product of vectors 1 and 2 to get normal vector of polygon face
    struct coordinate faceNormal;
    faceNormal.x = (AB.y * BC.z) - (AB.z * BC.y);
    faceNormal.y = (AB.z * BC.x) - (AB.x * BC.z);
    faceNormal.z = (AB.x * BC.y) - (AB.y * BC.x);

    return faceNormal;
}

/*
 * generateFaceNormalList(coordinateArray vertices, faceArray faces)
 * Input:
 * - vertices = A list of all the vertices of the polygons making up the polygon mesh for the sphere
 * - faces = A list of the faces of the polygon mesh as lists of the index of each vertex of the polygon in the vertices list
 * Output: Returns an array containing the normal vector for each polygon of the sphere polygon mesh
 * Description: Takes 3 points from a face by getting the index for the vertex from the face list and searching
 * for it in the vertices list. Sends these points to another function which uses these points to calculate the
 * surface normal for the polygon.
*/

coordinateArray generateFaceNormalList(coordinateArray vertices, faceArray faces)
{
    // Array to hold the calculated face normals
    coordinateArray normals;
    initCoordinateArray(&normals, 10);

    // Loop to get the vertex indices from faces, then to get the vertex coordinates matching each index from vertices
    for (int i = 0; i < faces.used; i++)
    {
        struct coordinate tempNormal = calculateFaceNormal(
            vertices.array[faces.array[i].array[0]],
            vertices.array[faces.array[i].array[1]],
            vertices.array[faces.array[i].array[2]]
        );

        tempNormal.distanceFromCamera = faces.array[i].distanceFromCamera;

        insertCoordinateArray(&normals, tempNormal);
    }

    return normals;
}

/*
 * calculateVectorBetweenTwoPoints(struct coordinate headPoint, struct coordinate tailPoint)
 * Input:
 * - headPoint = The coordinate at the front of the vector
 * - tailPoint = The coordinate at the back of the vector
 * Output: The components of the vector heading from tailPoint towards headPoint stored in a coordinate struct
 * Description: Uses two points to calculate the components of the vector heading from tailPoint towards
 * headPoint.
*/

struct coordinate calculateVectorBetweenTwoPoints(struct coordinate headPoint, struct coordinate tailPoint)
{
    struct coordinate vector;
    vector.x = headPoint.x = tailPoint.x;
    vector.y = headPoint.y = tailPoint.y;
    vector.z = headPoint.z = tailPoint.z;

    return vector;
}

/*
 * calculateVectorMagnitude(struct coordinate vector)
 * Input:
 * - vector = The vector stored as a coordinate
 * Output: The magnitude of the vector
 * Description: Uses the vector components to calculate its magnitude.
*/

float calculateVectorMagnitude(struct coordinate vector)
{
    return sqrt(pow(vector.x, 2.0) + pow(vector.y, 2.0) + pow(vector.z, 2.0));
}

/*
 * calculateVectorMagnitude(struct coordinate vector)
 * Input:
 * - vector = The vector stored as a coordinate
 * Output: The magnitude of the vector
 * Description: Uses the vector components to calculate its magnitude.
*/

float multiplyVectors(struct coordinate vector1, struct coordinate vector2)
{
    return (vector1.x * vector2.x) + (vector1.y * vector2.y) + (vector1.z * vector2.z);
}

/*
 * findAngleBetweenVectors(struct coordinate vector1, struct coordinate vector2)
 * Input:
 * - vector1 = The first vector
 * - vector2 = The second vector
 * Output: The angle between the vectors
 * Description: Multiplies the two vectors together, divides this value by their magnitudes
 * multiplied together and takes the inverse cosine of this final value to calculate the
 * angle between the two vectors.
*/

float findAngleBetweenVectors(struct coordinate vector1, struct coordinate vector2)
{
    float magnitudeOfV1 = calculateVectorMagnitude(vector1);
    float magnitudeOfV2 = calculateVectorMagnitude(vector2);

    float dotOfV1V2 = multiplyVectors(vector1, vector2);

    return acos(dotOfV1V2 / (magnitudeOfV1 * magnitudeOfV2));
}

/*
 * generateAmbientLight()
 * Input: N/A
 * Output: The amount of ambient light to be applied to a polygon as a coordinate representing the RGB values
 * Description: Multiplies the red, green and blue ambient light intensities by the corresponding material
 * coefficients of the polygon to get the intensity of each colour of light that will be used to render it.
*/

struct coordinate generateAmbientLight()
{
    struct coordinate rgb;
    rgb.x = AMBIENT_INTENSITY_R * SPHERE_AMBIENT_COEFFICIENT_R;
    rgb.y = AMBIENT_INTENSITY_G * SPHERE_AMBIENT_COEFFICIENT_G;
    rgb.z = AMBIENT_INTENSITY_B * SPHERE_AMBIENT_COEFFICIENT_B;

    return rgb;
}

/*
 * generateDiffuseLight(struct coordinate polygon, struct coordinate lightSource, struct coordinate polygonVector)
 * Input:
 * - polygon = The center point of the polygon
 * - lightSource = The coordinates of the light source
 * - polygonVector = The normal vector of the polygon
 * Output: The amount of diffuse light to be applied to a polygon as a coordinate representing the RGB values
 * Description: Calculates the vector from the light source to the center point of the polygon,
 * calculates the angle between this vector and the normal vector of the polygon, and uses these
 * values as well as the light intensities and material coefficients to calculate the intensity of
 * each colour of light to be rendered for the polygon.
*/

struct coordinate generateDiffuseLight(struct coordinate polygon, struct coordinate lightSource, struct coordinate polygonVector)
{
    struct coordinate lightSourceVector = calculateVectorBetweenTwoPoints(polygon, lightSource);

    float angleBetweenPVAndLSV = findAngleBetweenVectors(polygonVector, lightSourceVector);

    struct coordinate rgb;

    if (cos(angleBetweenPVAndLSV) < 0)
    {
        rgb.x = 0.0;
        rgb.y = 0.0;
        rgb.z = 0.0;
    }
    else
    {
        rgb.x = INTENSITY_R * SPHERE_DIFFUSE_COEFFICIENT_R * cos(angleBetweenPVAndLSV);
        rgb.y = INTENSITY_G * SPHERE_DIFFUSE_COEFFICIENT_G * cos(angleBetweenPVAndLSV);
        rgb.z = INTENSITY_B * SPHERE_DIFFUSE_COEFFICIENT_B * cos(angleBetweenPVAndLSV);
    }

    return rgb;
}

/*
 * generateSpecularLight(struct coordinate polygon, struct coordinate lightSource, struct coordinate camera)
 * Input:
 * - polygon = The center point of the polygon
 * - lightSource = The coordinates of the light source
 * - camera = The coordinates of the camera
 * Output: The amount of specular light to be applied to a polygon as a coordinate representing the RGB values
 * Description: Calculates the vector from the light source to the center point of the polygon,
 * calculates the vector from the center point of the polygon to the camera, calculates the angle
 * between these two vectors, and uses these values as well as the light intensities and material
 * coefficients to calculate the intensity of each colour of light to be rendered for the polygon.
*/

struct coordinate generateSpecularLight(struct coordinate polygon, struct coordinate lightSource, struct coordinate camera)
{
    struct coordinate lightSourceVector = calculateVectorBetweenTwoPoints(polygon, lightSource);
    struct coordinate cameraVector = calculateVectorBetweenTwoPoints(camera, polygon);

    float angleBetweenPVAndLSV = findAngleBetweenVectors(cameraVector, lightSourceVector);

    struct coordinate rgb;

    if (cos(angleBetweenPVAndLSV) < 0)
    {
        rgb.x = 0.0;
        rgb.y = 0.0;
        rgb.z = 0.0;
    }
    else
    {
        rgb.x = INTENSITY_R * SPHERE_SPECULAR_COEFFICIENT_R * cos(angleBetweenPVAndLSV);
        rgb.y = INTENSITY_B * SPHERE_SPECULAR_COEFFICIENT_G * cos(angleBetweenPVAndLSV);
        rgb.z = INTENSITY_B * SPHERE_SPECULAR_COEFFICIENT_B * cos(angleBetweenPVAndLSV);
    }

    return rgb;
}

/*
 * calculateTotalLightForPolygon(struct coordinate polygon, struct coordinate lightSource, struct coordinate polygonVector)
 * Input:
 * - polygon = The center point of the polygon
 * - lightSource = The coordinates of the light source
 * - polygonVector = The normal vector of the polygon
 * Output: The total amount of light to be applied to a polygon as a coordinate representing the RGB values
 * Description: Calculates the intensity values of all three types of light in red, green, and
 * blue and adds them together to return an RGB value representing the total amount of light
 * to be applied to the polygon for rendering.
*/

struct coordinate calculateTotalLightForPolygon(struct coordinate polygon, struct coordinate lightSource, struct coordinate polygonVector)
{
    struct coordinate ambient = generateAmbientLight();
    struct coordinate diffuse = generateDiffuseLight(polygon, lightSource, polygonVector);

    struct coordinate cameraCoordinates;
    cameraCoordinates.x = Ex;
    cameraCoordinates.y = Ey;
    cameraCoordinates.z = Ez;

    struct coordinate specular = generateSpecularLight(polygon, lightSource, cameraCoordinates);

    struct coordinate total;
    total.x = ambient.x + diffuse.x + specular.x;
    total.y = ambient.y + diffuse.y + specular.y;
    total.z = ambient.z + diffuse.z + specular.z;

    if (total.x > 1.0)
    {
        total.x = 1.0;
    }
    if (total.y > 1.0)
    {
        total.y = 1.0;
    }
    if (total.z > 1.0)
    {
        total.z = 1.0;
    }

    return total;
}

int main(int argc, char *argv[])
{  
    // Center point of the sphere
    point pt;

    pt.x = 0.0;
    pt.y = 0.0;
    pt.z = 0.0;

    struct coordinate lightSource;
    lightSource.x = Lx;
    lightSource.y = Ly;
    lightSource.z = Lz;

    // Fill the vertex, face, normal, and center point arrays to prepare to render the sphere

    coordinateArray a = generateVertexList(pt, SPHERE_RADIUS, LATITUDE_LINES, LONGITUDE_LINES);

    faceArray b = generateFaceList(LATITUDE_LINES, LONGITUDE_LINES);

    coordinateArray polygonCenterPoints = generatePolygonCenterPointList(a, b);

    coordinateArray c = generateFaceNormalList(a, b);

    // Sort the faces, normals, and center points in the order of furthest from the camera to closest so the
    // polygons are rendered in the correct order
    b = insertionSortFA(b);
    polygonCenterPoints = insertionSortCA(polygonCenterPoints);
    c = insertionSortCA(c);

    dmatrix_t E; /* The centre of projection for the camera */
    
    dmat_alloc(&E, 4, 1);
    
    E.m[1][1] = Ex;
    E.m[2][1] = Ey;
    E.m[3][1] = Ez;
    E.m[4][1] = 1.0;
    
    dmatrix_t G ; /* Point gazed at by camera */
    
    dmat_alloc(&G, 4, 1) ;
    
    G.m[1][1] = Gx;
    G.m[2][1] = Gy;
    G.m[3][1] = Gz;
    G.m[4][1] = 1.0;

    dmatrix_t C; /* The camera matrix */

    dmat_alloc(&C, 4, 4);
    C = *build_camera_matrix(&E, &G);

    // Setup the window and renderer
    SDL_Event event;
    SDL_Renderer *renderer;
    SDL_Window *window;
    int i;

    SDL_Init(SDL_INIT_VIDEO);
    SDL_CreateWindowAndRenderer(W, H, 0, &window, &renderer);
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);
    SDL_RenderClear(renderer);

    // Loop through the array of polygons
    for (i = 0; i < b.used; i++)
    {
        // Store the polygon vertices in matrices to be converted from world coordinates to screen coordinates
        dmatrix_t tempMatrixArray[b.array[i].used];

        for (int j = 0; j < b.array[i].used; j++)
        {
            dmatrix_t P;

            dmat_alloc(&P, 4, 1);

            P.m[1][1] = a.array[b.array[i].array[j]].x;
            P.m[2][1] = a.array[b.array[i].array[j]].y;
            P.m[3][1] = a.array[b.array[i].array[j]].z;
            P.m[4][1] = 1.0;

            P = *perspective_projection(dmat_mult(&C, &P));

            tempMatrixArray[j] = P;
        }
        
        //Calculate the light intensities used for each polygon and fill them with the provided filling algorithm
        struct coordinate totalLight = calculateTotalLightForPolygon(polygonCenterPoints.array[i], lightSource, c.array[i]);

        XFillConvexPolygon(totalLight, tempMatrixArray, renderer, b.array[i].used);
    }

    //Clean up resources before exiting
    SDL_RenderPresent(renderer);

    while (1)
    {
        if (SDL_PollEvent(&event) && event.type == SDL_QUIT)
        {
            break;
        }
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    freeCoordinateArray(&a);
    freeCoordinateArray(&c);
    freeFaceArray(&b);
}
