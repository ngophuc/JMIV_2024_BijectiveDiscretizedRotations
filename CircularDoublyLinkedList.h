#ifndef CircularDoublyLinkedList
#define CircularDoublyLinkedList

#include <vector>
#include "Tree.h"

using namespace std;

typedef struct {
  int idNode;
  vector<int> triplet;
  int norme;
  int alpha;
  int omega;
  nodeTree* p_segment;
} nodeData;

typedef struct node{
  nodeData data;
  struct node *next;
  struct node *prev;
} node;

node* getNodeFromId(node* head, int idNode);

node* createNode (int id, vector<int> t);

node* insertEnd(node** head, int id, vector<int> t);

node* insertAfter(node** head, int idNode, int id, vector<int> t);

node* updatePSegment(node* aNode, nodeTree* aNodeTree);

void deleteNode(node** head, int idNode);
void deleteNode(node** head, node* del_node);

void displayList(node* head);
void displayListBackward(node* head);
string displayNodeData(node* aNode);
void displayTree(node* head);

#endif // CircularDoublyLinkedList
