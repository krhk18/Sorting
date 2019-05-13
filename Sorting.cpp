/*
This program takes a integers and sort them, using different sorting algorithms:
- Insertion sort
- Shell sort
- Selection sort
- Merge sort
- Heap sort
- Quick sort

Kristina Håkansson, BTH 2019-04
*/

#include "Functions.h"
#include <iostream>
#include <vector>

int main()
{
	int data[] = { 13, 54, 16, 73, 5, 33, 21, 11, 75, 84 };
	int n = 10;	// (sizeof(data) / (sizeof(data[0]));
	
	std::cout << "Objects before sorting:\n";
	print(data, n);

	insertionSort(data, n);
	std::cout << "InsertionSort:\n";
	print(data, n);

	shellSort(data, n);
	std::cout << "ShellSort:\n";
	print(data, n);

	selectionSort(data, n);
	std::cout << "SelectionSort:\n";
	print(data, n);

	mergeSort(data, n);
	std::cout << "MergeSort:\n";
	print(data, n);

	heapSort(data, n);
	std::cout << "HeapSort:\n";
	print(data, n);

	quickSort(data, n);
	std::cout << "QuickSort:\n";
	print(data, n);

	std::cin.get();

	return 0;
}