#ifndef DoublyLinkedList
#define DoublyLinkedList

#include "Tree.h"

using namespace std;

typedef struct {
  int idNode;
  nodeTree* p_segment;
} nodeData;

typedef struct node{
  nodeData data;
  struct node *next;
  struct node *prev;
} node;

node* createNode (int id, nodeTree* aNodeTree);

node* insertAfter(node* prev_node, int id, nodeTree* aNodeTree);
node* insertEnd(node** head, int id, nodeTree* aNodeTree);
node* insertFront(node** head, int id, nodeTree* aNodeTree);

void deleteNode(node** head, node* del_node);

void displayList(node* head);

string displayNodeData(node* aNode);

node* getNodeFromId(node* head, int idNode);

void displayTree(node* head);

node* updatePSegment(node* aNode, nodeTree* aNodeTree);
#endif // DoublyLinkedList
