#include "linked_list_manager.h"

void print_list(node *head) {
    node *ptr = head;
    printf("[ ");
    while(ptr != NULL) {
        printf("(%d,%d,%lf) ",ptr->row,ptr->col,ptr->val);
        ptr = ptr->next;
    }
    printf("]\n");
}

void free_list(node *head)  {
    struct node* tmp;
    while (head != NULL)  {
        tmp = head;
        head = head->next;
        free(tmp);
    }
}

void insert_first(int row, int col, double val, node **head) {
    node *new_node = (node*) malloc(sizeof(node));
    new_node->row = row;
    new_node->col = col;
    new_node->val = val;
    
    new_node->next = *head;
    *head = new_node;
}

