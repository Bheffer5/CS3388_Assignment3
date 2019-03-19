/*
 * CS3388 Assignment 3
 * Author: Ben Heffernan
 * Date of Creation: March 6th, 2019
 * Purpose: A header file I used to hold all of the data structures I created for this project.
 * The code is pretty self explanatory as it just contains some simple structs, dynamic array
 * implementations using structs and pointers and sorting algorithms.
 */

#include <stdio.h>
#include <stdlib.h>

struct coordinate
{
    float x;
    float y;
    float z;
    float distanceFromCamera;
};

typedef struct
{
    float x;
    float y;
    float z;
}  point;

// Array implementation to hold lists of integers to represent polygons by the indices of their vertices in the vertex array

typedef struct
{
    int *array;
    float distanceFromCamera;
    size_t used;
    size_t size;
} vertexArray;

/*
 * initVertexArray(vertexArray *a, size_t initialSize)
 * Input:
 * - vertexArray = Integer pointer used to simulate a dynamic array of integers
 * - size_t = Initial size of the array
 * Output: void
 * Description: Allocates memory to initialize the pointer, sets the size of the array to initialSize
 * and sets the number of elements in the array (used) to 0.
*/

void initVertexArray(vertexArray *a, size_t initialSize)
{
    a->array = (int *)malloc(initialSize * sizeof(int));
    a->used = 0;
    a->size = initialSize;
}

/*
 * insertVertexArray(vertexArray *a, int element)
 * Input:
 * - vertexArray = Integer pointer used to simulate a dynamic array of integers
 * - element = Integer to add to the array
 * Output: void
 * Description: Adds the element to the array and resizes it if it is full.
*/

void insertVertexArray(vertexArray *a, int element)
{
    if (a->used == a->size)
    {
        a->size *= 2;
        a->array = (int *)realloc(a->array, a->size * sizeof(int));
    }

    a->array[a->used++] = element;
}

/*
 * freeVertexArray(vertexArray *a)
 * Input:
 * - vertexArray = Integer pointer used to simulate a dynamic array of integers
 * Output: void
 * Description: Frees the memory used by the array and sets the pointer to NULL
*/

void freeVertexArray(vertexArray *a)
{
    free(a->array);
    a->array = NULL;
    a->used = a->size = 0;
}

// Array implementation to hold the vertex arrays that represent the polygons

typedef struct
{
    vertexArray *array;
    size_t used;
    size_t size;
} faceArray;

/*
 * initFaceArray(faceArray *a, size_t initialSize)
 * Input:
 * - faceArray = VertexArray pointer used to simulate a dynamic array of vertexArrays
 * - size_t = Initial size of the array
 * Output: void
 * Description: Allocates memory to initialize the pointer, sets the size of the array to initialSize
 * and sets the number of elements in the array (used) to 0.
*/

void initFaceArray(faceArray *a, size_t initialSize)
{
    a->array = (vertexArray *)malloc(initialSize * sizeof(vertexArray));
    a->used = 0;
    a->size = initialSize;
}

/*
 * insertFaceArray(faceArray *a, int element)
 * Input:
 * - faceArray = VertexArray pointer used to simulate a dynamic array of vertexArrays
 * - element = Integer to add to the array
 * Output: void
 * Description: Adds the element to the array and resizes it if it is full.
*/

void insertFaceArray(faceArray *a, vertexArray element)
{
    if (a->used == a->size)
    {
        a->size *= 2;
        a->array = (vertexArray *)realloc(a->array, a->size * sizeof(vertexArray));
    }

    a->array[a->used++] = element;
}

/*
 * freeFaceArray(faceArray *a)
 * Input:
 * - faceArray = VertexArray pointer used to simulate a dynamic array of vertexArrays
 * Output: void
 * Description: Frees the memory used by the array and sets the pointer to NULL
*/

void freeFaceArray(faceArray *a)
{
    free(a->array);
    a->array = NULL;
    a->used = a->size = 0;
}

// Array implementation to hold all of the vertices as coordinate structs for the polygons of the meshes

typedef struct
{
    struct coordinate *array;
    size_t used;
    size_t size;
} coordinateArray;

/*
 * initCoordinateArray(coordinateArray *a, size_t initialSize)
 * Input:
 * - coordinateArray = Coordinate struct pointer used to simulate a dynamic array of coordinate structs
 * - size_t = Initial size of the array
 * Output: void
 * Description: Allocates memory to initialize the pointer, sets the size of the array to initialSize
 * and sets the number of elements in the array (used) to 0.
*/

void initCoordinateArray(coordinateArray *a, size_t initialSize)
{
    a->array = (struct coordinate *)malloc(initialSize * sizeof(struct coordinate));
    a->used = 0;
    a->size = initialSize;
}

/*
 * insertCoordinateArray(coordinateArray *a, int element)
 * Input:
 * - coordinateArray = Coordinate struct pointer used to simulate a dynamic array of coordinate structs
 * - element = Integer to add to the array
 * Output: void
 * Description: Adds the element to the array and resizes it if it is full.
*/

void insertCoordinateArray(coordinateArray *a, struct coordinate element)
{
    if (a->used == a->size)
    {
        a->size *= 2;
        a->array = (struct coordinate *)realloc(a->array, a->size * sizeof(struct coordinate));
    }

    a->array[a->used++] = element;
}

/*
 * freeCoordinateArray(coordinateArray *a)
 * Input:
 * - coordinateArray = Coordinate struct pointer used to simulate a dynamic array of coordinate structs
 * Output: void
 * Description: Frees the memory used by the array and sets the pointer to NULL
*/

void freeCoordinateArray(coordinateArray *a)
{
    free(a->array);
    a->array = NULL;
    a->used = a->size = 0;
}

/*
 * insertionSortCA(coordinateArray arr)
 * Input:
 * - coordinateArray = Coordinate struct array to be sorted
 * Output: The sorted coordinate struct array
 * Description: A simple implementation of the insertion sort algorithm to sort a coordinateArray
 * from furthest distance of the coordinate from the camera to the closest.
*/

coordinateArray insertionSortCA(coordinateArray arr)
{
    float max;

    for (int i = 0; i < arr.used; i++)
    {
        max = arr.array[i].distanceFromCamera;

        for (int j = i + 1; j < arr.used; j++)
        {
            if (arr.array[j].distanceFromCamera > max)
            {
                max = arr.array[j].distanceFromCamera;
                struct coordinate temp = arr.array[i];
                arr.array[i] = arr.array[j];
                arr.array[j] = temp;
            }
        }
    }

    return arr;
}

/*
 * insertionSortFA(faceArray arr)
 * Input:
 * - coordinateArray = VertexArray array to be sorted
 * Output: The sorted vertexArray array
 * Description: A simple implementation of the insertion sort algorithm to sort a faceArray
 * from furthest distance of the face from the camera to the closest.
*/

faceArray insertionSortFA(faceArray arr)
{
    float max;

    for (int i = 0; i < arr.used; i++)
    {
        max = arr.array[i].distanceFromCamera;

        for (int j = i + 1; j < arr.used; j++)
        {
            if (arr.array[j].distanceFromCamera > max)
            {
                max = arr.array[j].distanceFromCamera;
                vertexArray temp = arr.array[i];
                arr.array[i] = arr.array[j];
                arr.array[j] = temp;
            }
        }
    }

    return arr;
}
