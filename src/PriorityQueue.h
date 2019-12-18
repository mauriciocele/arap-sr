#if !defined(___PRIORITY_QUEUE___H___)
#define ___PRIORITY_QUEUE___H___

#include <vector>

/*
 * Implementation **based** in:
 *
 * Sedgewick Robert, "Priority Queues", Algorithms (pp.127-143), 1984
 * Addison-Wesley Series in Computer Science, ISBN O-201-06672-6
 *
 * The Heap data structure is like a binary tree arranged inside an array.
 * The root of the tree is at the zero position, the left child is at the 
 * k * 2 + 1 position and the right in the k * 2 + 2. On the other side
 * the parent of a child is in the ( k - 1 ) / 2 position.
 * The high priority value is always located at the root of the heap. We always insert
 * a new value at the bottom of the heap and compare it with its father, if it
 * is greater than its father then swap the values and repeat the operation again.
 * Similarly, we always delete the root of the heap, and place the bottom value as the
 * new root. Then we compare the root with its childs, if the childs are greater than
 * the root then we swap the value of the bigger child and repeat the operation again.
 */
template < typename T >
class PriorityQueue
{
private:
	//priority values labeled 0...N
	std::vector<T>		a;
	//index of priorities: a[ heap[k] ]
	std::vector<int>	heap;
	//inverse of index: inv[ heap[k] ] = heap[ inv[k] ]
	//	heap[x] = y
	//	inv[y] = x
	std::vector<int>	inv;
	//capacity of the heap
	int N;
	//number of actual elements on the heap
	int Count;
	//higher values have high priority (descending) or lower values have high priority (ascending)
	bool asc;

	/*
	 *	Return the parent of the k node. If the k node is the root then return -1.
	 */
	inline int		parent( int k ) { return ( k - 1 ) >> 1; }

	/*
	 *	Return the left child of the k node.
	 */
	inline int		child( int k ) { return ( k << 1 ) + 1; }

	/*
	 *	Move the k node down the tree.
	 *	Compare k with its childs, if the childs are greater than
	 *	k then swap the value of the bigger child and repeat the operation again.
	 */
	void downheap( int k )
	{
		int v = heap[k];
		while( k <= parent(Count-1) )
		{
			int j = child(k);

			if( j+1 < Count &&
				(a[ heap[j+1] ] > a[ heap[j] ] && !asc || 
				 a[ heap[j+1] ] < a[ heap[j] ] && asc) ) j++;
			if( a[v] >= a[ heap[j] ] && !asc ||
				a[v] <= a[ heap[j] ] && asc )
				break;
			heap[k] = heap[j];
			inv[heap[j]] = k;
			k = j;
		}
		heap[k] = v;
		inv[v] = k;
	}

	/*
	 *	Move the new inserted value up the tree.
	 *	Compare k with its father, if it is greater than its father 
	 *	then swap the values and repeat the operation again.
	 */
	void upheap( int k )
	{
		int v = heap[k];
		while( parent(k) >= 0 && 
			  (a[ heap[parent(k)] ] <= a[v] && !asc ||
			   a[ heap[parent(k)] ] >= a[v] && asc) )
		{
			heap[ k ] = heap[ parent(k) ];
			inv[ heap[ parent(k) ] ] = k;

			k = parent(k);
		}
		heap[ k ] = v;
		inv[v] = k;
	}

public:
	/*
	 *	Return the index of the root element.
	 */
	inline int		getTopIndex() { return heap[0]; }

	/*
	 *	Return the priority of the root element.
	 */
	inline T		getTopPriority() { return a[ heap[0] ]; }

	/*
	 *	Return the priority of the k element.
	 */
	inline T		getPriority( int k ) { return a[ k ]; }

	/*
	 *	Return true if the index is in the Priority Queue.
	 */
	inline bool		onHeap(int k) { return heap[ inv[ k ] ] == k; }

	/*
	 *	Return true if the Priority Queue is empty.
	 */
	inline bool		empty() { return Count == 0; }

	/*
	 *	Return the size of the Priority Queue
	 */
	inline int		size() { return Count; }

	/*
	 *	Return the capacity of the Priority Queue.
	 */
	inline int		capacity() { return N; }

	/*
	 *	N: total capacity of the queue
	 *	asc: true if lower values have high priority, false if higher values have high priority
	 */
	PriorityQueue( int N, bool asc )
	{
		this->N = N;
		this->asc = asc;
		Count = 0;

		heap.reserve( N );
		inv.reserve(N);
		a.reserve(N);

		for( int k = 0 ; k < N ; ++k )
		{
			heap.push_back( k );
			inv.push_back( k );
			a.push_back( 0 );
		}
	}

	/*
	 *	Construct the Priority Queue based on a predefined list of labeled priorities.
	 */
	void construct( const std::vector<T> &a )
	{
		for( int k = 0 ; k < N ; ++k )
			this->a[k] = a[k];

		Count = N;
		for( int k = parent(Count-1) ; k >= 0 ; --k )
			downheap( k );
	}

	/*
	 *	Insert a new value to the Priority Queue.
	 */
	void insert( int k, T priority )
	{
		a[ k ] = priority;
		heap[ Count ] = k;
		inv[ k ] = Count;
		Count++;
		upheap( Count-1 );
	}

	/*
	 *	Delete the root of the Priority Queue.
	 */
	void pop()
	{
		heap[0] = heap[Count-1];
		Count--;
		downheap( 0 );
	}

	/*
	 *	Remove an element from the Priority Queue.
	 */
	void remove( int k )
	{
		if( a[heap[Count-1]] < a[k] && !asc ||
			a[heap[Count-1]] > a[k] && asc )
		{
			heap[inv[k]] = heap[Count-1];
			Count--;
			downheap( inv[k] );
		}
		else
		{
			heap[inv[k]] = heap[Count-1];
			Count--;
			upheap( inv[k] );
		}
	}

	/*
	 *	Change the priority of an element.
	 */
	void changePriority( int k, T priority )
	{
		if( a[ k ] > priority && !asc ||
			a[ k ] < priority && asc )
		{
			a[ k ] = priority;
			downheap( inv[k] );
		}
		else
		{
			a[ k ] = priority;
			upheap( inv[k] );
		}
	}
};

///*
// * Implementation based in:
// *
// * Sedgewick Robert, "Priority Queues", Algorithms (pp.127-143), 1984
// * Addison-Wesley Series in Computer Science, ISBN O-201-06672-6
// *
// * The Heap data structure is like a binary tree arranged inside an array.
// * The root of the tree is at the zero position, the left child is at the 
// * k * 2 + 1 position and the right in the k * 2 + 2. On the other side
// * the parent of a child is in the ( k - 1 ) / 2 position.
// * The high priority value is always located at the root of the heap. We always insert
// * a new value at the bottom of the heap and compare it with its father, if it
// * is greater than its father then swap the values and repeat the operation again.
// * Similarly, we always delete the root of the heap, and place the bottom value as the
// * new root. Then we compare the root with its childs, if the childs are greater than
// * the root then we swap the value of the bigger child and repeat the operation again.
// */
//class PriorityQueue
//{
//private:
//	//priority values labeled 0...N
//	std::vector<int>	a;
//	//index of priorities: a[ heap[k] ]
//	std::vector<int>	heap;
//	//inverse of index: inv[ heap[k] ] = heap[ inv[k] ]
//	//	heap[x] = y
//	//	inv[y] = x
//	std::vector<int>	inv;
//	//capacity of the heap
//	int N;
//	//number of actual elements on the heap
//	int Count;
//
//	/*
//	 *	Return the parent of the k node. If the k node is the root then return -1.
//	 */
//	inline int		parent( int k ) { return ( k - 1 ) >> 1; }
//
//	/*
//	 *	Return the left child of the k node.
//	 */
//	inline int		child( int k ) { return ( k << 1 ) + 1; }
//
//	/*
//	 *	Move the k node down the tree.
//	 *	Compare k with its childs, if the childs are greater than
//	 *	k then swap the value of the bigger child and repeat the operation again.
//	 */
//	void downheap( int k )
//	{
//		int v = heap[k];
//		while( k <= parent(Count-1) )
//		{
//			int j = child(k);
//
//			if( j+1 < Count && a[ heap[j+1] ] > a[ heap[j] ] ) j++;
//			if( a[v] >= a[ heap[j] ] )
//				break;
//			heap[k] = heap[j];
//			inv[heap[j]] = k;
//			k = j;
//		}
//		heap[k] = v;
//		inv[v] = k;
//	}
//
//	/*
//	 *	Move the new inserted value up the tree.
//	 *	Compare k with its father, if it is greater than its father 
//	 *	then swap the values and repeat the operation again.
//	 */
//	void upheap( int k )
//	{
//		int v = heap[k];
//		while( parent(k) >= 0 && a[ heap[parent(k)] ] <= a[v] )
//		{
//			heap[ k ] = heap[ parent(k) ];
//			inv[ heap[ parent(k) ] ] = k;
//
//			k = parent(k);
//		}
//		heap[ k ] = v;
//		inv[v] = k;
//	}
//
//public:
//	/*
//	 *	Return the index of the root element.
//	 */
//	inline int		getTopIndex() { return heap[0]; }
//
//	/*
//	 *	Return the priority of the root element.
//	 */
//	inline int		getTopPriority() { return a[ heap[0] ]; }
//
//	/*
//	 *	Return the priority of the k element.
//	 */
//	inline int		getPriority( int k ) { return a[ k ]; }
//
//	/*
//	 *	Return true if the index is in the Priority Queue.
//	 */
//	inline bool		onHeap(int k) { return heap[ inv[ k ] ] == k; }
//
//	/*
//	 *	Return true if the Priority Queue is empty.
//	 */
//	inline bool		empty() { return Count == 0; }
//
//	/*
//	 *	Return the size of the Priority Queue
//	 */
//	inline int		size() { return Count; }
//
//	/*
//	 *	Return the capacity of the Priority Queue.
//	 */
//	inline int		capacity() { return N; }
//
//	PriorityQueue( int N )
//	{
//		this->N = N;
//		Count = 0;
//
//		heap.reserve( N );
//		inv.reserve(N);
//		a.reserve(N);
//
//		for( int k = 0 ; k < N ; ++k )
//		{
//			heap.push_back( k );
//			inv.push_back( k );
//			a.push_back( 0 );
//		}
//	}
//
//	/*
//	 *	Construct the Priority Queue based on a predefined list of labeled priorities.
//	 */
//	void construct( const std::vector<int> &a )
//	{
//		for( int k = 0 ; k < N ; ++k )
//			this->a[k] = a[k];
//
//		Count = N;
//		for( int k = parent(Count-1) ; k >= 0 ; --k )
//			downheap( k );
//	}
//
//	/*
//	 *	Insert a new value to the Priority Queue.
//	 */
//	void insert( int k, int priority )
//	{
//		a[ k ] = priority;
//		heap[ Count ] = k;
//		inv[ k ] = Count;
//		Count++;
//		upheap( Count-1 );
//	}
//
//	/*
//	 *	Delete the root of the Priority Queue.
//	 */
//	void pop()
//	{
//		heap[0] = heap[Count-1];
//		Count--;
//		downheap( 0 );
//	}
//
//	/*
//	 *	Remove an element from the Priority Queue.
//	 */
//	void remove( int k )
//	{
//		if( a[heap[Count-1]] < a[k] )
//		{
//			heap[inv[k]] = heap[Count-1];
//			Count--;
//			downheap( inv[k] );
//		}
//		else
//		{
//			heap[inv[k]] = heap[Count-1];
//			Count--;
//			upheap( inv[k] );
//		}
//	}
//
//	/*
//	 *	Change the priority of an element.
//	 */
//	void changePriority( int k, int priority )
//	{
//		if( a[ k ] > priority )
//		{
//			a[ k ] = priority;
//			downheap( inv[k] );
//		}
//		else
//		{
//			a[ k ] = priority;
//			upheap( inv[k] );
//		}
//	}
//};

#endif