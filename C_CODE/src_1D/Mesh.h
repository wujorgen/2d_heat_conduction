#ifndef MESH_H
#define MESH_H

#ifndef NODE_H
#include "Node.h"
#endif

struct Mesh
{
    struct Node Nodes[]; // can we use calloc this in main?
};

#endif