public class BinarySearch_ {

    // if no duplicates present  -->  a, b, c
    //if duplicates present      --> a(will return any of found position), d, e

    // a
    // for array with duplicate elements, it might return any of the found index.
    public static int binarySearch(int[] nums, int target) {
        int n = nums.length;
        int low = 0, high = n - 1;
        while (low <= high) {
            int mid = (low + high) / 2;
            if (nums[mid] == target) return mid;
            else if (target > nums[mid]) low = mid + 1;
            else high = mid - 1;
        }
        return -1;
    }

    //b  --> searchInsert or Ceil
    // if target exists return its index else return the index of element greater than target
    public int ceil(int[] nums, int target) {
        int l=0,h=nums.length;
        while(l<h){
            int mid = l+(h-l)/2;
            if(nums[mid]>=target)h=mid;
            else l=mid+1;
        }
        return l;
    }
    // c
    // if target exists return its index else return the index of element smaller than target
    public int floor(int[] nums, int target) {
        int l=0,h=nums.length;
        while(l<h){
            int mid = l+(h-l)/2;
            if(nums[mid]>target)h=mid;
            else l=mid+1;
        }
        return l-1;
    }

    // d  -> lower bound
    // find the exact target's smallest index in the array of duplicate elements
    // same code as ceil (** not same code as floor)

    // e  -> upper bound
    // find the exact target's largest index in the array of duplicate elements
    // same code as floor (** not same code as ceil)

}
