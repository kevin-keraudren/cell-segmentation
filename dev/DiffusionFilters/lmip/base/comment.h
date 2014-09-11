/*--------------------------------------------------------------------------*/
/** \file comment.h

   \brief Basic handling of comments for data files. 

   Maintainer: Stephan Didas <didas@mia.uni-saarland.de>                    */
/*--------------------------------------------------------------------------*/

#ifndef __COMMENT_H__
#define __COMMENT_H__

#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/*--------------------------------------------------------------------------*/

/** \struct comment
    \brief Represents a comment as a simply connected list of strings. 
 */
struct comment {
    char*  line;           /**< one line of text          */
    struct comment *next;  /**< pointer to the next entry */
};

/*--------------------------------------------------------------------------*/

/** \brief Allocates storage for a comment. 

    The pointers line and next are initialised to NULL.
    This is the smallest list except the empty one. 
 */
void comment_alloc(struct comment **c)    /**< pointer to the comment */
{
    *c = (struct comment*) malloc( sizeof( struct comment ) );

    (*c)->line = NULL;
    (*c)->next = NULL;

    return;
}

/*--------------------------------------------------------------------------*/

/** \brief Disallocates storage of a comment.

    Since comments are simply connected lists, this function disallocates 
    all elements in the list after the argument c and including c.
 */
void comment_disalloc(struct comment *c)    /**< the comment */
{
    struct comment *next;

    if(c == NULL)
    {
	return;
    }

    while(c->next != NULL)
    {
	next = c->next;

	if(c->line != NULL)
	{
	    free(c->line);
	}
	free(c);

	c = next;
    }

    if(c->line != NULL)
    {
	free(c->line);
    }
    free(c);

    return;
}

/*--------------------------------------------------------------------------*/

/** \brief Adds one line to a comment.

    This function has a variable number of arguments such that it can be 
    used in the same way as the printf function family. 

    For example, it is possible to write 
       comment_add(c, "number of iterations %ld\n", it);
    to involve variable values in the comment. 

    \note The argument c has already to be allocated.
 */
void comment_add(struct comment *c,    /**< the comment */
		 const char* format,   /**< format string for the new line */
		 ...)                  /**< further arguments (see above) */
{
    char* line;
    struct comment *new_comment;
    va_list args;
    va_start(args, format);

    if(c == NULL)
    {
	printf("comment_add: comment is NULL.\n");
	exit(1);
    }

    line = (char*) malloc(256 * sizeof(char));

    vsnprintf(line, 255, format, args); 

    if((c->line) == NULL) 
    {
	c->line = line;
    } else {
	/* go to the end of the list */
	while(c->next != NULL)
	{
	    c = c->next;
	}

	/* create a new node and add it */
	comment_alloc(&new_comment);
	new_comment->line = line;
	c->next = new_comment;
    }

    va_end(args);

    return;
} /* comment_add */

/*--------------------------------------------------------------------------*/

/** \brief Adds a separator line to a comment.
 */
void comment_add_sep(struct comment* c)
{
    comment_add(c, "# -------------------------------------------------"
		"--------------------\n");

    return;
}

/*--------------------------------------------------------------------------*/

/** \brief Adds one line with the current time to a comment.
 */
void comment_add_time(struct comment* c)
{
    time_t now;

    now = time(0);
    comment_add(c, "# time: %s", ctime(&now) );

    return;
}

/*--------------------------------------------------------------------------*/

/** \brief Appends one comment to another one.

    Remember that comments are lists: This function appends one of these 
    lists to another one. 
    This is necessary to append the older comments to the new ones in 
    an application, such that in an output file the newest comments are on
    top.
    Be careful: Since the lists are appended, disallocating com1 will also 
    delete com2. 
    If you want to keep com2, use commend_append_copy instead.
 */
void comment_append(
    struct comment *com1,  /**< the first comment, the second one will be 
			        appended here
			    */
    struct comment *com2)  /**< the second comment */
{
    if(com1 == NULL) {
	com1 = com2;
    }

    while(com1->next != NULL)
    {
	com1 = com1->next;
    }

    com1->next = com2;

    return;
} /* comment_append */

/*--------------------------------------------------------------------------*/

/** \brief Copies one comment into another one.

    \note The destination comment has to be allocated before.
    If the source comment is NULL, the destination one is disallocated 
    also. 
 */
void comment_copy(
    struct comment *source, /**< the source comment      */
    struct comment *dest)   /**< the destination comment */
{
    long size;
    char *line;

    if(source == NULL) 
    {
	printf("source comment is NULL.\n");
	comment_disalloc(dest);
    }

    while(source != NULL && source->line != NULL)
    {
	size = strlen(source->line);
	line = (char*)malloc((size + 1) * sizeof(char));
	strcpy(line, source->line);
	dest->line = line;

	if(source->next != NULL) {
	    dest->next = (struct comment*)malloc(sizeof(struct comment));
	    dest   = dest->next;
	} else {
	    dest->next = NULL;
	}
	source = source->next;
    }

    return;
} /* comment_copy */

/*--------------------------------------------------------------------------*/

/** \brief Appends a copy of one comment to another one.

    This is useful is you want to use one comment in multiple files
    with slight modifications.
    For example, one can read the comment of one input file and 
    append it to the comments of all output files. 
 */
void comment_append_copy(
    struct comment *com1,  /**< the first comment, a copy of the second 
			        one will be appended here
			    */
    struct comment *com2)  /**< the second comment */
{
    struct comment *help;

    if(com1 == NULL) {
	printf("comment 1 is zero.\n");
	comment_alloc(&com1);
	comment_copy(com2, com1);
	return;
    }

    while(com1->next != NULL)
    {
	com1 = com1->next;
    }
    comment_alloc(&help);
    comment_copy(com2, help);

    com1->next = help;

    return;
} /* comment_append_copy */

/*--------------------------------------------------------------------------*/

/** \brief Prints a comment out on screen.
   
    \note The argument c has already to be allocated.
 */
void comment_print(struct comment *c)
{
    if(c == NULL) {
	printf("comment_print: NULL comment!\n");
	exit(1);
    }

    if(c->line == NULL) {
	printf("comment_print: comment is empty\n");
    } else {
	while(c != NULL)
	{
	    if(c->line != NULL) 
	    {
		printf("%s", c->line);
	    }
	    c = c->next;
	}
    }

    return;
}

/*--------------------------------------------------------------------------*/

/** \brief Prints a comment out in a stream.
   
    \note The argument c has already to be allocated.
 */
void comment_fprint(FILE *stream, struct comment *c)
{
    if(c == NULL) {
	printf("comment_print: NULL comment!\n");
	return;
    }

    if(c->line == NULL) {
	printf("comment_print: comment is empty\n");
    } else {
	while(c != NULL)
	{
	    if(c->line != NULL)
	    {
		fprintf(stream, "%s", c->line);
	    }
	    c = c->next;
	}
    }

    return;
}

/*--------------------------------------------------------------------------*/

/** \brief Replaces one line in a comment.

    This can be useful if a program writes multiple output files which have
    most of the parameters in common. One only has to replace the lines 
    for which the information has changed.

    \note Starts to count the lines with zero.
 */
void comment_replace_line(
    struct comment* c,    /**< comment                        */
    long   number,        /**< number of the line to replace  */
    const char* format,   /**< format string for the new line */
    ...)                  /**< further arguments (see above)  */
{
    char* line;
    long  i = 0;

    va_list args;
    va_start(args, format);

    if(c == NULL)
    {
	printf("comment_replace_line: comment is NULL.\n");
	exit(1);
    }

    line = (char*) malloc(256 * sizeof(char));

    vsnprintf(line, 255, format, args); 

    while(c->next != NULL && i < number)
    {
	c = c->next;
	i = i+1;
    }
    
    if(i != number) 
    {
	printf("comment_replace_line: should replace line no. %ld while\n",
	       number);
	printf("comment_replace_line: comment only has %ld lines\n", i);

	exit(1);
    }

    /* if there is an old line, free the memory */
    if(c->line != NULL) 
    {
	free(c->line);
    } 

    /* hang the new line in */
    c->line = line;

    va_end(args);
     
    return;
} /* comment_replace_line */

/*--------------------------------------------------------------------------*/

/** \brief Replaces one line of a comment with the current time.
 */
void comment_replace_time(
    struct comment* c, /**< the comment */
    long number)       /**< number of the line */
{
    time_t now;

    now = time(0);
    comment_replace_line(c, number, "# time: %s", ctime(&now) );

    return;
}

/*--------------------------------------------------------------------------*/

#endif /* __COMMENT_H__ */

