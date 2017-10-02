///////////////////////////////////////////////////////////////////////////////
//                   ALL STUDENTS COMPLETE THESE SECTIONS
// Main Class File:  TestSort.java
// File:             TestSort.java, ComparisonSort.java, 
//					 SortObject.java
// Semester:         CS367 Summer 2016
//
// Author:           Jason Choe choe2@wisc.edu
// CS Login:         choe
// Lecturer's Name:  Amanda Strominger
// Lab Section:      001
//
//////////////////// PAIR PROGRAMMERS COMPLETE THIS SECTION ////////////////////
//
// Pair Partner:     Jason Leduc
// Email:            jlleduc@wisc.edu
// CS Login:         leduc
// Lecturer's Name:  Amanda Strominger
// Lab Section:      001
//
////////////////////////////////////////////////////////////////////////////////

/**
 * This class implements six different comparison sorts (and an optional
 * seventh sort for extra credit):
 * 
 * 1. selection sort
 * 2. insertion sort
 * 3. merge sort
 * 4. quick sort
 * 5. heap sort
 * 6. selection2 sort
 * (extra credit) 
 * 7. insertion2 sort
 * 
 * It also has a method that runs all the sorts on the same input array and
 * prints out statistics.
 * 
 * There is also a helper method "swap" to exchange the contents of two indices
 * in an array
 * 
 * Bugs: None known
 *
 * @author Jason Choe, Jason Leduc
 */

public class ComparisonSort {
	// dataMoves field counts reassignments in the array
	private static int dataMoves;

    /**
     * Sorts the given array using the selection sort algorithm discussed in 
     * lecture
     * 
     * @param E[]A   		 the array to sort
     */
    public static <E extends Comparable<E>> void selectionSort(E[] A) {
    	if(A == null || A.length == 0 ) return;
    	for(int i = 0 ; i < A.length ; i++){
    		int minIndex = i;
    		for(int j = i+1; j < A.length  ; j++){
    			if(A[j].compareTo(A[minIndex]) < 0){
    				minIndex = j;
    			}
    		}
    		// if next leftmost spot doesn't already contain the new minimum 
    		// value, swap the next leftmost element in the unsorted (right) 
    		// part of the array with the minimum value found in the unsorted 
    		// part.
    		if(i != minIndex){
    			swap(A, i, minIndex);
    		}
    	}
    }

    /**
     * Sorts the given array using the insertion sort algorithm. Note: after
     * this method finishes the array is in sorted order.
     * 
     * @param [E]A    the array to sort
     */
    public static <E extends Comparable<E>> void insertionSort(E[] A) {
    	if(A == null || A.length == 0 ) return;
    	for(int i = 1 ; i < A.length ; i++){
    		// A[0] to A[i-1] is sorted
    		E temp = A[i];
    		dataMoves++;
    		int j;
    		for(j = i - 1; j >= 0 && A[j].compareTo(temp) > 0 ; j--){
    			A[j+1] = A[j]; // move A[j] over to the right
    			dataMoves++;
    		}
    		A[j+1] = temp;
    		dataMoves++;
    	}
    }

    /**
     * Sorts the given array using the merge sort algorithm. Note: after this
     * method finishes the array is in sorted order.
     * 
     * @param E[]A    the array to sort
     */
    
	public static <E extends Comparable<E>> void mergeSort(E[] A) {
    	if(A == null || A.length == 0 ) return;
	    mergeAux(A, 0, A.length - 1); // call the aux. function to do	
	    							  // all the work
	}
	
	/**
	 * mergeSort helper method
	 * @param E[]A		the array to sort
	 * @param low		index of left end of array being sorted
	 * @param high		index of right end of array being sorted
	 */ 
	private static <E extends Comparable<E>> void mergeAux(E[] A, int low, 
			int high) {
		
	    // base case
	    if (low == high) return;
	 
	    // recursive case
	    
	 // Step 1: Find the middle of the array (conceptually, divide it in half)
	    int mid = (low + high) / 2;
	     
	 // Steps 2 and 3: Sort the 2 halves of A
	    mergeAux(A, low, mid);
	    mergeAux(A, mid+1, high);
	 
	 // Step 4: Merge sorted halves into an auxiliary array
	    E[] tmp = (E[])(new Comparable[high-low+1]);
	    int left = low;    // index into left half
	    int right = mid+1; // index into right half
	    int pos = 0;       // index into tmp
	     
	    while ((left <= mid) && (right <= high)) {
	    	// choose the smaller of the two values "pointed to" by left, right
	    	// copy that value into tmp[pos]
	    	// increment either left or right as appropriate
	    	// increment pos

    	    if (A[left].compareTo(A[right]) <= 0) {
    	        tmp[pos] = A[left];
    	        dataMoves++;
    	        left++;
    	    }
    	    else {
    	        tmp[pos] = A[right];
    	        dataMoves++;
    	        right++;
    	    }
    	    pos++;
    	}
	     // when one of the two sorted halves has "run out" of values, but
		 // there are still some in the other half; copy all the remaining 
		 // values to tmp
		 // Note: only 1 of the next 2 loops will actually execute
		 while (left <= mid) {
		     tmp[pos] = A[left];  
		     dataMoves++;
		     left++;
		     pos++;
		 }
		 while (right <= high) {
		     tmp[pos] = A[right]; 
		     dataMoves++;
		     right++;
		     pos++;
		 }       
	    // all values are in tmp; copy them back into A
		 for(int i = 0 ; i < tmp.length ; i++){
			 A[i+low] = tmp[i];
			 dataMoves++;
		 }
	}

    /**
     * Sorts the given array using the quick sort algorithm, using the median of
     * the first, last, and middle values in each segment of the array as the
     * pivot value. Note: after this method finishes the array is in sorted
     * order.
     * 
     * @param E[]A   the array to sort
     */
	public static <E extends Comparable<E>> void quickSort(E[] A) {
    	if(A == null || A.length == 0 ) return;
	    quickSort(A, 0, A.length-1);
	}
	 
	/**
	 * quickSort helper for recursive sorting of partitions
	 * @param E[]A		the array to sort
	 * @param left		index of left end
	 * @param right		index of right end
	 */
    private static <E extends Comparable<E>> void quickSort(E[] A, int left, 
    			int right) {
    	int index = partition(A, left, right);
        if (left < index - 1)
        	quickSort(A, left, index - 1);
        if (index < right)
            quickSort(A, index, right);
    }
	 
    /**
	 * Partition method of quick sort which sorts based on a pivot determined by
	 * setPivot
	 * 
	 * @param E[]A		the array to sort
	 * @param left	    the index of the low end of array to quick sort
	 * @param right     the index of the high end of array to quick sort
	 * @return position of where the swapped pivot was moved to
	 */
    private static <E extends Comparable<E>> int partition(E[] A, int left, 
    		int right){
        int i = left+1, j = right-2;
        setPivot(A, left, right);
        E pivot = A[right-1];
        dataMoves++;
        // start by comparing 2nd element and 3rd from end with pivot	
        // and work your way in until the left pointer finds a value greater
        // than the pivot and the right pointer finds a value less than the 
        // pivot, then we swap them, and continue this process until the 
        // pointers cross
        while (i <= j) {
            while (A[i].compareTo(pivot) < 0)
                  i++;
            while (A[j].compareTo(pivot) > 0)
                  j--;
            if (i <= j) {
            	swap(A, i, j);
                  i++;
                  j--;
            }
        }
        return i;
    }
    
    /**
     * setPivot compares the left and right values and the middle value to 
     * put the smallest at the leftmost position in the array, the largest in
     * the rightmost position, and the median one just left of the rightmost, 
     * where
     * it will act as the pivot
     * @param E[]A		the array to sort
     * @param left		index of left end
     * @param right		index of right end
     */
    private static <E extends Comparable<E>> void setPivot(E[] A, int left, 
    		int right){
    	int middle = (left + right) / 2;
		if(A[left].compareTo(A[middle]) > 0)
			swap(A, left, middle);
		if(A[middle].compareTo(A[right]) > 0){
			swap(A, middle, right);
			if(A[left].compareTo(A[middle]) > 0)
				swap(A, left, middle);
		}
		swap(A, middle, right-1);
    }
    
    /**
     * Sorts the given array using the heap sort algorithm outlined below. Note:
     * after this method finishes the array is in sorted order.
     * <p>
     * The heap sort algorithm is:
     * </p>
     * 
     * <pre>
     * for each i from 1 to the end of the array
     *     insert A[i] into the heap (contained in A[0]...A[i-1])
     *     
     * for each i from the end of the array up to 1
     *     remove the max element from the heap and put it in A[i]
     * </pre>
     * 
     * @param E[]A    the array to sort
     */
    public static <E extends Comparable<E>> void heapSort(E[] A) {
    	if(A == null || A.length == 0 ) return;
    	E temp; // for temporarily storing items when swapping positions
    	E[] heapArray = (E[])(new Comparable[A.length + 1]);

    	// insert all elements of A into heapArray
    	// heap starts at index 1 so make the array one longer than A
    	for(int i = 1 ; i < heapArray.length ; i++){
    		heapArray[i] = A[i-1];
    		dataMoves++;
    		// need to modify the index in the while loop below so make copy
    		int index = i;
    		// heap keeps minimum element at the front (index 1)
    		while(index > 1 && heapArray[index].compareTo(heapArray[index/2]) <
    				0){
    			temp = heapArray[index];
    			heapArray[index] = heapArray[index/2];

    			heapArray[index/2] = temp;
    			dataMoves += 3;
    			index = index/2;
    		}
    	}

    	// remove all elements from heapArray in ascending order
    	// keep track of where the last element in the heap is in the array
    	// and rearrange the heap to keep the minimum element at the front
    	int numItems = heapArray.length - 1;
    	for(int i = 0 ; i < A.length ; i++){

    		A[i] = heapArray[1];
    		heapArray[1] = heapArray[numItems];
    		heapArray[numItems] = null;
    		dataMoves += 3;
    		numItems--;
    		int index = 1;
    		int swapIndex = 0;
    		while(index * 2 <= numItems){
    			if(heapArray[index * 2].compareTo(heapArray[index]) < 0){
    				swapIndex = index * 2;
    			}
    			if(index * 2 + 1 <= numItems)
    				if((heapArray[index * 2 + 1].compareTo(heapArray[index]) 
    						< 0) && (heapArray[index * 2 + 1].compareTo
    								(heapArray[index * 2]) < 0))
    					swapIndex = index * 2 + 1;
    			if(swapIndex < index * 2){
    				// set index so while loop terminates
    				index = numItems;
    			}
    			else {
    				temp = heapArray[index];
    				heapArray[index] = heapArray[swapIndex];
    				heapArray[swapIndex] = temp;
    				dataMoves += 3;
    				index = swapIndex;
    			}
    		}
    	}
    }
    		
     /**
     * Sorts the given array using the selection2 sort algorithm outlined
     * below. Note: after this method finishes the array is in sorted order.
     * 
     * The selection2 sort is a bi-directional selection sort that sorts
     * the array from the two ends towards the center. The selection2 sort
     * algorithm is:
     * 
     * begin = 0, end = A.length-1
     * 
     * At the beginning of every iteration of this loop, we know that the 
     * elements in A are in their final sorted positions from A[0] to 
     * A[begin-1] and from A[end+1] to the end of A.  That means that 
     * A[begin] to A[end] are still to be sorted.
     * do use the MinMax algorithm (described below) to find the minimum and 
     * maximum values between A[begin] and A[end]
     *     
     *     swap the maximum value and A[end]
     *     swap the minimum value and A[begin]
     *     
     *     ++begin, --end
     * until the middle of the array is reached
     * 
     * The MinMax algorithm allows you to find the minimum and maximum of N
     * elements in 3N/2 comparisons (instead of 2N comparisons). The way to do
     * this is to keep the current min and max; then
     * 
     * take two more elements and compare them against each other
     * compare the current max and the larger of the two elements
     * compare the current min and the smaller of the two elements
     * 
     * @param E[]A    the array to sort
     */
    public static <E extends Comparable<E>> void selection2Sort(E[] A) {
    	if(A == null || A.length == 0 ) return;
    	int newMaxIndex, newMinIndex, begin, end;
    	E temp;
    	// gradually move min to the front and max to the end
    	for(begin = 0, end = A.length - 1 ; begin < A.length/2 ; begin++, 
    			end--) {
    		newMinIndex = begin;
    		newMaxIndex = end; 
			// find newMin and newMax to replace current begin and end values
    		for(int j = begin , k = end ; j <= A.length/ 2  ; j++, k--){
    			// compare the two new ends with each other
    			if(A[j].compareTo(A[k]) < 0){
    				// compare lesser end with current minimum to possibly 
    				// replace
    				if(A[j].compareTo(A[newMinIndex]) < 0 )
    					newMinIndex = j;
    				// compare greater end with current maximum
    				if(A[k].compareTo(A[newMaxIndex]) > 0)
    					newMaxIndex = k;
    			} else {
    				if(A[k].compareTo(A[newMinIndex]) < 0 )
    					newMinIndex = k;
    				if(A[j].compareTo(A[newMaxIndex]) > 0)
    					newMaxIndex = j;
    			}
    		}

    		// handle special cases where max or min start out at the 
    		// beginning or end.  This adds code but saves some datamoves 
    		// over the alternative method of comparing the ends and swapping 
    		//them to put the smallest first and the largest at the end
    		if(begin == newMaxIndex && end == newMinIndex)
    			swap(A, begin, end);
    		else if(begin == newMinIndex && end == newMaxIndex){
    			continue;
    		}
    		else if(begin == newMaxIndex){
    			swap(A, begin, end);
    			swap(A, begin, newMinIndex);
    		} else if(end == newMinIndex){
    			swap(A, begin, end);
    			swap(A, end, newMaxIndex);
    		} else {
    			swap(A, begin, newMinIndex);
    			swap(A, end, newMaxIndex);
    		}
    	}
    }

    
    /**
     * Extra Credit: Sorts the given array using the insertion2 sort 
     * algorithm outlined below.  Note: after this method finishes the array 
     * is in sorted order.
     * 
     * The insertion2 sort is a bi-directional insertion sort that sorts the 
     * array from the center out towards the ends.  The insertion2 sort 
     * algorithm is:
     *
     * precondition: A has an even length
     * left = element immediately to the left of the center of A
     * right = element immediately to the right of the center of A
     * if A[left] > A[right]
     *     swap A[left] and A[right]
     * left--, right++ 
     *  
     * At the beginning of every iteration of this loop, we know that the 
     * elements in A from A[left+1] to A[right-1] are in relative sorted order.
     * do
     *     if (A[left] > A[right])
     *         swap A[left] and A[right]
     *  
     * starting with with A[right] and moving to the left, use insertion sort 
     * algorithm to insert the element at A[right] into the correct location 
     * between A[left+1] and A[right-1]
     *     
     *  starting with A[left] and moving to the right, use the insertion sort 
     *  algorithm to insert the element at A[left] into the correct location 
     *  between A[left+1] and A[right-1]
     *  
     *  left--, right++
     * until left has gone off the left edge of A and right has gone off the 
     * right edge of A
     * 
     * This sorting algorithm described above only works on arrays of even 
     * length.  If the array passed in as a parameter is not even, the method 
     * throws an IllegalArgumentException
     *
     * @param  E[]A 	the array to sort
     * @throws IllegalArgumentException if the length or A is not even
     */
    public static <E extends Comparable<E>> void insertion2Sort(E[] A) { 
    	if(A == null || A.length == 0 ) return;
    	// the array must have an even number of elements
    	if(A.length%2 != 0)
    		throw new 
    		IllegalArgumentException("Must provide an even-sized array");
    	E temp; // temporary storage for edge values
    	int left = A.length/2 -1; // left center
    	int right = A.length/2; // right center
    	
    	// if left is greater than right, swap them
    	if(A[left].compareTo(A[right]) > 0)
    		swap(A, left, right);
    	// moves edges out
    	left--;
    	right++;

    	// loop with edges moving out after insertions
    	while(left >= 0){
    		// Checking right and left for swap
    		if(A[left].compareTo(A[right]) > 0)
    			swap(A, left, right);
    		int i = right;
    		// store the right to insert to proper place later
    		temp = A[right];
    		dataMoves++;
    		// Inserting from right
    		while(temp.compareTo(A[i-1]) < 0 && i > left){
    			A[i] = A[i-1];
    			dataMoves++;
    			i--;
    		}
    		// if right value was not greatest, move it to new location
    		if(A[right] != temp){
    			A[i] = temp;
    			dataMoves++;
    		}
    		i = left;
    		temp = A[left];
    		dataMoves++;
    		// inserting from left
    		while(temp.compareTo(A[i+1]) > 0 && i < right){
    			A[i] = A[i+1];
    			dataMoves++;
    			i++;
    		}
    		if(A[left] != temp){
    			A[i] = temp;
    			dataMoves++;
    		}
    		// move edges out
    		left--;
    		right++;
    	}
    }

    /**
     * Internal helper for printing rows of the output table.
     * 
     * @param sort          name of the sorting algorithm
     * @param compares      number of comparisons performed during sort
     * @param moves         number of data moves performed during sort
     * @param milliseconds  time taken to sort, in milliseconds
     */
    private static void printStatistics(String sort, int compares, int moves,
                                        long milliseconds) {
        System.out.format("%-23s%,15d%,15d%,15d\n", sort, compares, moves, 
                          milliseconds);
    }
    
    // used to look at arrays to see how sorting is going
    private static void print(Comparable[] arr){
    	int size;
    	if(arr.length < 15)
    		size = arr.length;
    	else
    		size = 15;
    	for(int i = 0 ; i < size ; i++){
    		if(arr[i] != null)
    			System.out.print(arr[i].toString() + ", ");
    	}
    	System.out.println();
    }
    
    /**
     *  simple method used to swap objects in two indices in an array
     * @param E[]A	the array to sort
     * @param i		index of item that will be switched
     * @param j		index of another item that will be switched
     */
    private static <E extends Comparable<E>> void swap (E[] A, int i, int j) {
        E temp = A[i];
        A[i] = A[j];
        A[j] = temp;
        dataMoves += 3 ;
    }

    /**
     * Sorts the given array using the six (seven with the extra credit)
     * different sorting algorithms and prints out statistics. The sorts 
     * performed are:
     * <ul>
     * <li>selection sort</li>
     * <li>insertion sort</li>
     * <li>merge sort</li>
     * <li>quick sort</li>
     * <li>heap sort</li>
     * <li>selection2 sort</li>
     * <li>(extra credit) insertion2 sort</li>
     * </ul>
     * <p>
     * The statistics displayed for each sort are: number of comparisons, 
     * number of data moves, and time (in milliseconds).
     * </p>
     * <p>
     * Note: each sort is given the same array (i.e., in the original order) 
     * and the input array A is not changed by this method.
     * </p>
     * 
     * @param A  the array to sort
     */
    static public void runAllSorts(SortObject[] A) {
        System.out.format("%-23s%15s%15s%15s\n", "algorithm", "data compares", 
                          "data moves", "milliseconds");
        System.out.format("%-23s%15s%15s%15s\n", "---------", "-------------", 
                          "----------", "------------");
        SortObject[] B;	//to store clone of A
        
        // selection sort
        B = A.clone(); // clone A so we can reuse it for other sorts
        SortObject.resetCompares();
        dataMoves = 0;
        // set initial timer
        long initialTime = System.currentTimeMillis();
        // data moves are returned from sort method
        selectionSort(B);
        // get time after sort completed
        long finalTime = System.currentTimeMillis();
        // check comparison counter
        int dataCompares = SortObject.getCompares();
        // check that the data are sorted
        for(int i = 0 ; i < B.length - 1 ; i++)
        	if(B[i].compareTo(B[i+1]) > 0){
        		System.out.println("NOT SORTED!!!");
        		System.exit(1);
        	}
        printStatistics("selection", dataCompares, dataMoves, finalTime -
        		initialTime);
        
        // insertion sort
        B = A.clone();
        SortObject.resetCompares();
        dataMoves = 0;
        initialTime = System.currentTimeMillis();
        insertionSort(B);
        finalTime = System.currentTimeMillis();
        dataCompares = SortObject.getCompares();
        // check that it's sorted
        for(int i = 0 ; i < B.length - 1 ; i++){
        	if(B[i].compareTo(B[i+1]) > 0){
        		System.out.println("NOT SORTED!!!");
        		System.exit(1);
        	}
        }
        printStatistics("insertion", dataCompares, dataMoves, finalTime -
        		initialTime);
        
        // merge sort
        B = A.clone();
        SortObject.resetCompares();
        dataMoves = 0;
        initialTime = System.currentTimeMillis();
        mergeSort(B);
        finalTime = System.currentTimeMillis();
        dataCompares = SortObject.getCompares();
        // check that it's sorted
        for(int i = 0 ; i < B.length - 1 ; i++){
        	if(B[i].compareTo(B[i+1]) > 0){
        		System.out.println("NOT SORTED!!!");
        		System.exit(1);
        	}
        }
        printStatistics("merge", dataCompares, dataMoves, finalTime -
        		initialTime);
        
        // QuickSort
        B = A.clone();
        SortObject.resetCompares();
        dataMoves = 0;
        initialTime = System.currentTimeMillis();
        quickSort(B);
        finalTime = System.currentTimeMillis();
        dataCompares = SortObject.getCompares();
        // check that it's sorted
        for(int i = 1 ; i < B.length ; i++){
        	if(B[i-1].compareTo(B[i]) > 0){
        		System.out.println("NOT SORTED!!!");
        		System.exit(1);
        	}
        }
        printStatistics("quick", dataCompares, dataMoves, finalTime -
        		initialTime);
        
        // HeapSort
        B = A.clone();
        SortObject.resetCompares();
        dataMoves = 0;
        initialTime = System.currentTimeMillis();
        heapSort(B);
        finalTime = System.currentTimeMillis();
        dataCompares = SortObject.getCompares();
        // check that it's sorted
        for(int i = 1 ; i < B.length ; i++){
        	if(B[i-1].compareTo(B[i]) > 0){
        		System.out.println("NOT SORTED!!!");
        		System.exit(1);
        	}
        }
        printStatistics("heap", dataCompares, dataMoves, finalTime -
        		initialTime);
        
        // Selection2Sort
        B = A.clone();
        SortObject.resetCompares();
        dataMoves = 0;
        initialTime = System.currentTimeMillis();
        selection2Sort(B);
        finalTime = System.currentTimeMillis();
        dataCompares = SortObject.getCompares();
        // check that it's sorted
        for(int i = 1 ; i < B.length ; i++){
        	if(B[i-1].compareTo(B[i]) > 0){
        		System.out.println("NOT SORTED!!!");
        		System.exit(1);
        	}
        }
        printStatistics("selection2", dataCompares, dataMoves, finalTime -
        		initialTime);
        
        
        // insertion2Sort
        B = A.clone();
        SortObject.resetCompares();
        dataMoves = 0;
        initialTime = System.currentTimeMillis();
        insertion2Sort(B);
        finalTime = System.currentTimeMillis();
        dataCompares = SortObject.getCompares();
        // check that it's sorted
        for(int i = 1 ; i < B.length ; i++){
        	if(B[i-1].compareTo(B[i]) > 0){
        		System.out.println("NOT SORTED!!!");
        		System.exit(1);
        	}
        }
        printStatistics("insertion2", dataCompares, dataMoves, finalTime -
        		initialTime);
    }
}

