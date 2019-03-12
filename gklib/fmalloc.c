/*
  * allocate N bytes of memory on-the-fly
*/

#include <stdlib.h>

#if defined(ESV) || defined(OSX)
#else
#include <malloc.h>
#endif

#if defined(HPUX) || defined(VMS)
#   define FMALLOC fmalloc		/* VMS, HPUX need no trailing _  */
#   define FFREE ffree
#   define void int
#else
#   define FMALLOC fmalloc_	/* but ESV, SGI, LINUX, ALPHA, OSX do */
#   define FFREE ffree_
#endif

#ifdef ALPHA
typedef long address_type;
#else
typedef int address_type;
#endif


/*
  *  fmalloc -- allocate memory dynamically.
  *  Called from FORTRAN: adr = fmalloc(size). The address can then be
  *  passed to a subprogram like this: call sub(%val(adr), size),
  *  and fmalloc should be declared as integer (*8 on ALPHA !)
  */
address_type *FMALLOC(size)
      int *size;
{
   address_type *ptr;
   size_t n;
/*  char *malloc(); */

/*   n = (*size) * sizeof(int); */

   n = (*size);
   ptr = (address_type *)malloc( n );
   return ptr;
}

/*
  * ffree -- deallocate dynamic memory that was previously allocated
  * by a call to fmalloc. Called from FORTRAN: call ffree (adr)
  */
void FFREE(ptr)
      address_type *ptr;
{
     free( (void *)*ptr);

}
