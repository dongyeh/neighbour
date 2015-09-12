#ifndef GLOBAL_H
#define GLOBAL_H

typedef int NeighbourPair[2];

#define PAIR_FACTOR 100
#define FLT_EPSILON 1e-8
#define TREE_INCREMENT 3
#define NEGATIVE_INFINITY -1e20

#define tnt tree[now]
#define tlson tree[tnt.lson]
#define trson tree[tnt.rson]
#define tleaf tree[tnt.toLeaf]

#endif // GLOBAL_H
