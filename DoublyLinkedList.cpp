#include <iostream>

#include "DoublyLinkedList.h"

node* createNode (int id, nodeTree* aNodeTree) {
  node* aNode = new node();
  // assign data to the node
  aNode->data.idNode = id;
  aNode->data.p_segment = aNodeTree;
  // assign prev pointer to the node
  aNode->prev = NULL;
  // assign next pointer to the node
  aNode->next = NULL;
  
  return aNode;
}

// insert a node after a specific node
node* insertAfter(node* prev_node, int id, nodeTree* aNodeTree) {

  // check if previous node is NULL
  if (prev_node == NULL) {
    cout << "previous node cannot be NULL";
    return;
  }
  
  // allocate memory for newNode
  node* newNode = createNode(id, aNodeTree);
  
  // set next of newNode to next of prev node
  newNode->next = prev_node->next;
  
  // set next of prev node to newNode
  prev_node->next = newNode;
  
  // set prev of newNode to the previous node
  newNode->prev = prev_node;
  
  // set prev of newNode's next to newNode
  if (newNode->next != NULL)
    newNode->next->prev = newNode;
  
  return newNode;
}

// insert a newNode at the end of the list
node* insertEnd(node** head, int id, nodeTree* aNodeTree) {
  // allocate memory for newNode
  node* newNode = createNode(id, aNodeTree);
  
  // assign NULL to next of newNode
  newNode->next = NULL;
  
  // store the head node temporarily (for later use)
  node* temp = *head;
  
  // if the linked list is empty, make the newNode as head node
  if (*head == NULL) {
    newNode->prev = NULL;
    *head = newNode;
    return newNode;
  }
  
  // if the linked list is not empty, traverse to the end of the linked list
  while (temp->next != NULL)
    temp = temp->next;
  
  // now, the last node of the linked list is temp
  
  // point the next of the last node (temp) to newNode.
  temp->next = newNode;
  
  // assign prev of newNode to temp
  newNode->prev = temp;
  
  return newNode;
}

// insert node at the front
node* insertFront(node** head, int id, nodeTree* aNodeTree) {
  // allocate memory for newNode
  node* newNode = createNode(id, aNodeTree);
  
  // point next of newNode to the first node of the doubly linked list
  newNode->next = (*head);
  
  // point prev to NULL
  newNode->prev = NULL;
  
  // point previous of the first node (now first node is the second node) to newNode
  if ((*head) != NULL)
    (*head)->prev = newNode;
  
  // head points to newNode
  (*head) = newNode;
  
  return newNode;
}

// delete a node from the doubly linked list
void deleteNode(node** head, node* del_node) {
  // if head or del is null, deletion is not possible
  if (*head == NULL || del_node == NULL)
    return;

  // if del_node is the head node, point the head pointer to the next of del_node
  if (*head == del_node)
    *head = del_node->next;

  // if del_node is not at the last node, point the prev of node next to del_node to the previous of del_node
  if (del_node->next != NULL)
    del_node->next->prev = del_node->prev;

  // if del_node is not the first node, point the next of the previous node to the next node of del_node
  if (del_node->prev != NULL)
    del_node->prev->next = del_node->next;

  // free the memory of del_node
  free(del_node);
}

// print the doubly linked list
string displayVector(std::vector<int> v) {
  return "(" + std::to_string(v[0]) + ", "+ std::to_string(v[1]) + ", "+ std::to_string(v[2]) + ")" ;
}

string displayNodeData(node* aNode) {
  string s = "[" + std::to_string(aNode->data.idNode) + ", " + displayVector(aNode->data.p_segment->tripletLeft) + "--" + displayVector(aNode->data.p_segment->tripletRight) + "]";
  return s;
}

void displayList(node* head) {
  node* last = NULL;

  while (head != NULL) {
    //cout << node->data << "->";
    cout << displayNodeData(head) << "->";
    last = head;
    head = head->next;
  }
  if (head == NULL)
    cout << "NULL\n";
}

void displayTree(node* head)
{
  node* temp = head;
  
  //printf("\nTraversal in forward direction \n");
  while (head!=NULL && temp->next != head) {
    nodeTree* n = temp->data.p_segment;
    vector<int> left = n->tripletLeft;
    vector<int> right = n->tripletRight;
    cout << displayVector(left) << "->" << displayVector(right) <<endl;
    temp = temp->next;
  }
  if(temp!=NULL) {
    nodeTree* n = temp->data.p_segment;
    vector<int> left = n->tripletLeft;
    vector<int> right = n->tripletRight;
    cout << displayVector(left) << "->" << displayVector(right) <<endl;
  }
  else
    cout << "Begin --> End: Head" << endl;
}

node* updatePSegment(node* aNode, nodeTree* aNodeTree){
  aNode->data.p_segment = aNodeTree;
  return aNode;
}
