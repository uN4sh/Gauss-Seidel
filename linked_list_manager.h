#ifndef __linked_list_h
#define __linked_list_h

#include <stdio.h>
#include <stdlib.h>

typedef struct node {
   int row;
   int col;
   double val;
   struct node *next;
} node;

void print_list(node *head);

void free_list(node *head);

void insert_first(int row, int col, double val, node **head);

#endif