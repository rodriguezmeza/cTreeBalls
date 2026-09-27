/* Typed introsort for the exact lexicographic neighbor order. Inlining the
 * comparison avoids qsort's indirect callback and generic record swaps.
 * Median-of-three partitioning has a heapsort depth limit; recursion always
 * follows the smaller partition, so adversarial input cannot exhaust stack.
 * Ordinal is a unique final key, preserving qsort's deterministic total order. */
static inline int lya_neighbor_less(const lya_neighbor *a,const lya_neighbor *b)
{
    if(a->forest_id!=b->forest_id) return a->forest_id<b->forest_id;
    if(a->leg_bin!=b->leg_bin) return a->leg_bin<b->leg_bin;
    if(a->radius!=b->radius) return a->radius<b->radius;
    return a->ordinal<b->ordinal;
}

static void lya_neighbor_sift(lya_neighbor *a,size_t root,size_t count)
{
    lya_neighbor value=a[root];
    while(root<count/2) {
        size_t child=2*root+1;
        if(child+1<count && lya_neighbor_less(a+child,a+child+1)) child++;
        if(!lya_neighbor_less(&value,a+child)) break;
        a[root]=a[child];root=child;
    }
    a[root]=value;
}

static void lya_neighbor_heap(lya_neighbor *a,size_t count)
{
    for(size_t i=count/2;i>0;i--) lya_neighbor_sift(a,i-1,count);
    for(size_t i=count;i>1;i--) {
        lya_neighbor value=a[0];a[0]=a[i-1];a[i-1]=value;
        lya_neighbor_sift(a,0,i-1);
    }
}

static void lya_neighbor_intro(lya_neighbor *a,size_t count,unsigned depth)
{
    while(count>16) {
        if(!depth) {lya_neighbor_heap(a,count);return;}
        depth--;
        const lya_neighbor *x=a,*y=a+count/2,*z=a+count-1,*temp;
        if(lya_neighbor_less(y,x)){temp=x;x=y;y=temp;}
        if(lya_neighbor_less(z,y)){temp=y;y=z;z=temp;}
        if(lya_neighbor_less(y,x)){temp=x;x=y;y=temp;}
        lya_neighbor pivot=*y;
        size_t left=0,right=count-1;
        for(;;) {
            while(lya_neighbor_less(a+left,&pivot)) left++;
            while(lya_neighbor_less(&pivot,a+right)) right--;
            if(left>=right) break;
            lya_neighbor value=a[left];a[left]=a[right];a[right]=value;
            left++;right--;
        }
        /* The median pivot lies within the array; nonempty partitions are
         * guaranteed by the total order. Keep a defensive bounded fallback. */
        if(left==0 || left==count) {lya_neighbor_heap(a,count);return;}
        if(left<count-left) {
            lya_neighbor_intro(a,left,depth);a+=left;count-=left;
        } else {
            lya_neighbor_intro(a+left,count-left,depth);count=left;
        }
    }
    for(size_t i=1;i<count;i++) {
        lya_neighbor value=a[i];size_t j=i;
        while(j && lya_neighbor_less(&value,a+j-1)) {a[j]=a[j-1];j--;}
        a[j]=value;
    }
}

static void lya_sort_neighbors(lya_neighbor *a,size_t count)
{
    unsigned depth=0;for(size_t n=count;n>1;n/=2)depth+=2;
    lya_neighbor_intro(a,count,depth);
}
