#include <iostream>

#include "CircularDoublyLinkedList.h"

node* createNode (int id, vector<int> t/*vector<int> tf, const vector<int> tr*/) {
  node* aNode = new node();
  // assign data to the node
  aNode->data.idNode = id;
  aNode->data.triplet = t;
  aNode->data.norme = t[0]*t[0]+t[1]*t[1];
  aNode->prev = NULL;
  // assign next pointer to the node
  aNode->next = NULL;
  
  return aNode;
}

// Function to insert at the end
node* insertEnd(node** head, int id, vector<int> t)
{
  // If the list is empty, create a single node
  // circular and doubly list
  if (*head == NULL) {
    node* new_node = createNode(id, t);
    new_node->next = new_node->prev = new_node;
    *head = new_node;
    return NULL;
  }
  
  // If list is not empty
  
  /* Find last node */
  node* last = (*head)->prev;
  
  // Create Node dynamically
  node* new_node = createNode(id, t);
  // head is going to be next of new_node
  new_node->next = *head;
  
  // Make new node previous of head
  (*head)->prev = new_node;
  
  // Make last previous of new node
  new_node->prev = last;
  
  // Make new node next of old last
  last->next = new_node;
  
  return new_node;
}

// Function to insert new node after the node with idNode
node* insertAfter(node** head, int idNode, int id, vector<int> te)
{
  // Create Node dynamically
  node* new_node = createNode(id, te);
  
  // Find node having idNode and next node of it
  node* temp = *head;
  while (temp->data.idNode != idNode)
    temp = temp->next;
  node* next = temp->next;
  
  // insert new_node between temp and next.
  temp->next = new_node;
  new_node->prev = temp;
  new_node->next = next;
  next->prev = new_node;
  
  return new_node;
}

// Function to delete a given node from the list
void deleteNode(node** head, int idNode) {
    // If list is empty
    if (*head == NULL)
        return;
  
    // Find the required node
    // Declare two pointers and initialize them
    node *curr = *head, *prev_1 = NULL;
    while (curr->data.idNode != idNode) {
        // If node is not present in the list
        if (curr->next == *head) {
            printf("\nList doesn't have node with value = %d", idNode);
            return;
        }
  
        prev_1 = curr;
        curr = curr->next;
    }
  
    // Check if node is the only node in list
    if (curr->next == *head && prev_1 == NULL) {
        (*head) = NULL;
        free(curr);
        return;
    }
  
    // If list has more than one node,
    // check if it is the first node
    if (curr == *head) {
        // Move prev_1 to last node
        prev_1 = (*head)->prev;
  
        // Move head ahead
        *head = (*head)->next;
  
        // Adjust the pointers of prev_1 and head node
        prev_1->next = *head;
        (*head)->prev = prev_1;
        free(curr);
    }
  
    // check if it is the last node
    else if (curr->next == *head) {
        // Adjust the pointers of prev_1 and head node
        prev_1->next = *head;
        (*head)->prev = prev_1;
        free(curr);
    }
    else {
        // create new pointer, points to next of curr node
        node* temp = curr->next;
  
        // Adjust the pointers of prev_1 and temp node
        prev_1->next = temp;
        temp->prev = prev_1;
        free(curr);
    }
}

// delete a node from the doubly linked list
void deleteNode(node** head, node* del_node) {
  deleteNode(head, del_node->data.idNode);
}


// print the doubly linked list
string displayVector(std::vector<int> v) {
  return "(" + std::to_string(v[0]) + ", "+ std::to_string(v[1]) + ", "+ std::to_string(v[2]) + ")" ;
}

string displayNodeData(node* aNode) {
  string s = "[" + std::to_string(aNode->data.idNode) + ", " + displayVector(aNode->data.triplet) + "]";
  return s;
}


void displayList(node* head)
{
  node* temp = head;
  
  //printf("\nTraversal in forward direction \n");
  while (head!=NULL && temp->next != head) {
    //printf("%d ", temp->data);
    cout << displayNodeData(temp) << "->";
    temp = temp->next;
  }
  if(temp!=NULL)
    //printf("%d ", temp->data);
    cout << displayNodeData(temp) << "--> End: Head" << endl;
  else
    cout << "Begin --> End: Head" << endl;
}

void displayTree(node* head)
{
  node* temp = head;
  
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

void displayListBackward(node* head)
{
  node* temp = head;
  
  node* last = head->prev;
  temp = last;
  while (temp->prev != last) {
    cout << displayNodeData(temp) << "->";
    temp = temp->prev;
  }

  cout << displayNodeData(temp) << "--> End: Head" << endl;
}

node* getNodeFromId(node* head, int idNode){
  // Find node having idNode and next node of it
  node* temp = head;
  while (temp->data.idNode != idNode)
    temp = temp->next;
  //node* next = temp->next;//element after idNode
  return temp;
}

node* updatePSegment(node* aNode, nodeTree* aNodeTree){
  aNode->data.p_segment = aNodeTree;
  return aNode;
}

