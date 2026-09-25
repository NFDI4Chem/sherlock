#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "dyntext.h"

/* Global definitions */

#define EOL (-1)
#define SIZE_CUR_STRING 5
#define MAX_SIZE_FRAG_NAME 31

/* token identifiers */
enum token_value {
	AND=0, OR, NOT, FRAG, FRNM, LP, RP, ERR, END 
} ;

/* opcode mnemonics for the Logic Virtual Machine */ 
enum opcodes_mn {
	JZ=0, JNZ, FR, LBL, INV
} ;

/* opcode definition, for forward and backward chaining */
typedef struct Opcode {
/* pointer to next opcode */
	struct Opcode *next ;
/* pointer to previous opcode */
	struct Opcode *prev ;
/* pointer to nexp opcode for JZ an JNZ opcodes,
to be used when the test succeeds */
	struct Opcode *jump ;
/* opcode mnemonic, as defined by enum opcodes_mn */
	int mnemo ;
/* address, used for jump, label and fragment opcodes */
	int addr ;
/* value for fragment opcodes, caches fragment search */
	int value ;
/* adress of a value, points to the single value that is
concerned by a given fragment addr(ess) */
	int *pvalue ;
} opcode ;
typedef opcode *popcode ;
#define ONE ((popcode)1)

/* a lvm structure contains the compilation and execution
environment of a logical expression */
typedef struct Lvm {
/* the string that constains the logical expression
that will be parsed before evaluation.
Il must be initialized before anything is started. */
	char *parsed_string ;

/* cursor for parsed_string, initial value is -1.
Reading parsed_string increments this index, so that
the pointed character is always the current one. */
	int parsed_string_index ;

/* value of last parsed token id */
	int cur_token ;

/* stores current string for variable id and logical operators */
	char cur_string[MAX_SIZE_FRAG_NAME] ;

/* pointers to first and last opcodes */
	popcode root ;
	popcode last ;

/* number of independent variables */
	int fragc ;

/* the vector of variable (numeric) names */
	int *fragv ;

/* counts calls to get_token() */
#ifdef DEBUG
	int ngt ;
#endif
} lvm ;
typedef lvm *plvm ;

/* Evaluation function type definition */

typedef int (*lvm_eval)(int) ;

/* function prototypes */

plvm lvm_new(char *s) ;
int lvm_error(plvm pl) ;
void lvm_vars(plvm pl, int *pfragc, int **pfragv) ;
int lvm_run(plvm pl, lvm_eval f) ;
void lvm_clean(plvm pl) ;

static int gc(plvm pl) ;
static void pb(plvm pl) ;
static int eol(plvm pl) ;
static int get_token(plvm pl) ;
static void add_code(plvm pl, int op, int label) ;
static void print_code (popcode p) ;
static void chop_last(plvm pl, int op, int num) ;
static void print_codes(plvm pl) ;
static void free_codes(plvm pl) ;
static void expr(plvm pl, int level) ;
static void term(plvm pl, int level) ;
static void prim(plvm pl, int level) ;
static popcode jump_search(plvm pl) ;
static void jump_proc(plvm pl) ;
static popcode frag_search(plvm pl) ;
static void frag_proc(plvm pl) ;
static int frag_value(popcode p, lvm_eval f) ;

extern int	fldeff ; /* lsd.c */
extern int flmaxfrag ; /* lsd.c */
extern dyntext	*strtab ; /* lsd.c */
extern dynstring	*curstring ; /* lsd.c */
extern void	test_open_skel(int f, int p) ; /* sub.c */

/* Class variables */

/* logical operator string values, followed by token string values */
static char *ops[] = {
	"AND", "OR", "NOT", "FRAG", "LP", "RP", "ERR", "END"
} ;

/* printable version of opcodes */
static char *opcodes_txt[] = {"JZ", "JNZ", "FRAG", "LBL", "NOT"} ;

/*  (The Only) Class method */

plvm lvm_new(char *s)
/* compiles string s, returns pointer to lvm object, NULL otherwise */
{
	plvm pl ;

/* memory allocation */
	pl = (plvm)malloc(sizeof(lvm)) ;
	if(pl == NULL) {
		return NULL ;
	}
/* initialisation of lvm members */
	pl->parsed_string_index = -1 ;
	pl->cur_token = END ;
	pl->root = NULL ;
	pl->last = NULL ;
	pl->fragc = 0 ;
	pl->fragv = NULL ;
#ifdef DEBUG
	pl->ngt = 0 ;
#endif

	pl->parsed_string = s ;
	expr(pl, 1) ;
#ifdef DEBUG
print_codes(pl) ;
#endif
	if(pl->cur_token != ERR) {
		jump_proc(pl) ;
		frag_proc(pl) ;
	}
	return pl ;
}

/* Public object methods */

int lvm_error(plvm pl)
/* returns 1 if an error occured, 0 otherwise */
{
	if(pl == NULL) return 1 ;
	return pl->cur_token == ERR ;	
}

void lvm_vars(plvm pl, int *pfragc, int **pfragv)
/* exports the number of variable and their numerical names */
{
	if(pl == NULL) return ;
	*pfragc = pl->fragc ;
	*pfragv = pl->fragv ;
}

int lvm_run(plvm pl, lvm_eval f)
/* simulates the Logical Virtual Machine */
{
	int v ;
	popcode p ;

/* resets all value field to -1 to set all values to
the undefined state */
	for(p=pl->root ; p != NULL ; p=p->next) {
		p->value = -1 ;
	}

/* processes opcodes*/
	p = pl->root ;
	while(p != NULL) {
#ifdef DEBUG
printf("Running opcode ") ;
print_code(p) ;
#endif
		switch(p->mnemo) {
		case FR :
			v = frag_value(p, f) ;
#ifdef DEBUG
printf("Value is now %d\n", v) ;
#endif
			p = p->next ;
			break ;
		case INV :
			v = !v ;
#ifdef DEBUG
printf("Value is now %d\n", v) ;
#endif
			p = p->next ;
			break ;
		case LBL :
			p = p->next ;
			break ;
		case JNZ :
			p = v ? p->jump : p->next ;
			break ;
		case JZ :
			p = v ? p->next : p->jump ;
			break ;
		}
	}
	return v ;
}

void lvm_clean(plvm pl)
/* call free_codes() and autodestruction */
{
	if(pl->fragv != NULL) free((char *)pl->fragv) ;
	free_codes(pl) ;
	free((char *)pl) ;
}

/* Private object methods */

static int gc(plvm pl)
/* gets next character in parsed_string.
returns EOL at end of string. */ 
{
	char c ;

	pl->parsed_string_index++ ;
	if ((c = pl->parsed_string[pl->parsed_string_index]) == '\0') {
#ifdef DEBUG
printf("EOL reached.\n") ;
#endif
		return EOL ;
	} else {
#ifdef DEBUG
printf("reading character %c\n", c) ;
#endif
		return (int)c ;
	}
}

static void pb(plvm pl)
/* moves cursor one position back in parsed_string,
silently does nothing if the cursor is at the beginning
of parsed_string */
{
	if(pl->parsed_string_index >= 0) {
#ifdef DEBUG
printf("putting back character %c\n",
	pl->parsed_string[pl->parsed_string_index]) ;
#endif
		pl->parsed_string_index-- ;
	}
}

static int eol(plvm pl)
/* returns 1 if end_of_line was reached, 0 otherwise */
{
	if(pl->parsed_string_index < 0) return 0 ;
	return pl->parsed_string[pl->parsed_string_index] == '\0' ;
}

static int get_token(plvm pl)
/* returns id of next token and sets cur_token
to the same value */
{
	int ic ;
	char c ;
	char *p ;
	int i ;
	int v ;
#ifdef DEBUG
	pl->ngt++ ;
#endif
/* skips heading blank characters */
	do {
		if((ic = gc(pl)) == EOL) {
			return pl->cur_token = END ;
		}
	} while (isspace(c = (char)ic)) ;
/* c is now the first non space character */
/* token analysis */
	switch(c) {
/* left parenthesis */
	case '(': 
#ifdef DEBUG
printf("token: Left parenthesis\n") ;
#endif
		return pl->cur_token = LP ;
/* right parenthesis */
	case ')':
#ifdef DEBUG
printf("token: Right parenthesis\n") ;
#endif
		return pl->cur_token = RP ;
/* fragment */
	case 'F':
	case 'f':
		p = pl->cur_string ;
		i = 0 ;
/* a fragment id is at most 4 digits long */
		while(
			i<SIZE_CUR_STRING &&
			(ic=gc(pl)) != EOL && 
			isdigit(c=(char)ic)) {
				p[i++] = c ;
		}
		pb(pl) ;
		if(i == SIZE_CUR_STRING) {
			return pl->cur_token = ERR ;
		}
		p[i] = '\0' ;
#ifdef DEBUG
printf("Fragment %s\n", p) ;
#endif
		return pl->cur_token = FRAG ;
/* fragment name */
	case '`':
		i = 0 ;
		dynstring_reset(curstring) ;
/* a fragment name is at most 31 characters long */
		while(
			i<MAX_SIZE_FRAG_NAME &&
			(ic=gc(pl)) != EOL && 
			isgraph(c=(char)ic) &&
			c != '`' ) {
				dynstring_pushc(curstring, c) ;
				i++ ;
		}
		if(i == MAX_SIZE_FRAG_NAME) {
			return pl->cur_token = ERR ;
		}
		flmaxfrag++ ;
		fldeff++ ;
		
		v = dyntext_push(strtab, curstring, 0, NULL) ;
/* assignment of an unused number of fragment Fn-like to the current fragment name */
		dyntext_setindex(strtab, v, flmaxfrag) ;
		test_open_skel(flmaxfrag, v) ; /* errors 431 to 436 */
		
#ifdef DEBUG
printf("Fragment name: %s\n", curstring->string) ;
#endif
		return pl->cur_token = FRNM ;		
/* logical operator */
	default :
		pb(pl) ;
		p = pl->cur_string ;
		i = 0 ;
/* a logical operator is at most 3 characters long */	
		while(
			i<3 &&
			(ic=gc(pl)) != EOL &&
			isalpha(c=(char)ic)) {
				c = toupper(c) ;
/* changes all lowercase to uppercase for compiling */
				p[i++] = c ;
		}
		p[i] = '\0' ;
/* looking for one of 3 logical operators */
		i = 0 ;
		while(i<3 && strcmp(p, ops[i])) {
			i++ ;
		}
#ifdef DEBUG
printf("Logical operator %d: %s %s\n", i, p, ops[i]) ;
#endif
		return pl->cur_token = ((i==3) ? ERR : i) ;
	}
}

static void add_code(plvm pl, int op, int num)
/* adds an opcode to the list pointer by root
and that ends by last. */
{
	popcode p ;

	p = (popcode)malloc(sizeof(opcode)) ;
	if(p == NULL) {
		pl->cur_token = ERR ;
		return ;
	}
	p->next = NULL ;
	p->jump = ONE ;
	p->mnemo = op ;
	p->addr = num ;
	p->value = -1 ;
	p->pvalue = NULL ;
	if(pl->root == NULL) {
		pl->root = p ;
		p->prev = NULL ;
	} else {
		pl->last->next = p ;
		p->prev = pl->last ;
	}
	pl->last = p ;
#ifdef DEBUG
printf("add ") ;
print_code(p) ;
#endif
}

static void print_code (popcode p)
/* prints a single opcode. */
{
	int i ;

	printf("%p: ", (void *)p) ;
	printf("%s ", opcodes_txt[p->mnemo]) ;
	if((i = p->addr) >= 0) {
		printf("%d", i) ;
	}
	printf("\n") ;
}

static void chop_last(plvm pl, int op, int num)
{
	popcode p ;

	if(pl->last == NULL && pl->last->mnemo != op && pl->last->addr != num) {
/* should never happen ! */
		pl->cur_token = ERR ;
	}
#ifdef DEBUG
printf("rem ") ;
print_code(pl->last) ;
#endif
	p = pl->last->prev ;
	free((char *)pl->last) ;
	if(p == NULL) {
		pl->root = pl->last = NULL ;
	} else {
		p->next = NULL ;
		pl->last = p ;
	}
}

static void print_codes(plvm pl)
/* print all opcodes, starting from root. */
{
	popcode p ;

	if(pl->root == NULL) return ;
	for(p = pl->root ; p != NULL ; p = p->next) {
		print_code(p) ;
	}
}

static void free_codes(plvm pl)
/* frees memory for all opcodes starting from root */
{
	popcode p ;

	while(pl->last != NULL) {
		p = pl->last->prev ;
		free(pl->last) ;
		pl->last = p ;
	}
}

static void expr(plvm pl, int level)
/* compiles an expression. An expression is a term, eventually ORed
with other terms. */
{
	int full_expr = 0 ;

	if(level == 1) {
		get_token(pl) ;
#ifdef DEBUG
printf("\t%d: %s\n", pl->ngt, ops[pl->cur_token]) ;
#endif
	}
	if(pl->cur_token == ERR) return ;
	term(pl, level+1) ;
	add_code(pl, JNZ, level) ;
	for(;;) {
		switch(pl->cur_token) {
		case OR :
			full_expr = 1 ;
			get_token(pl) ;
#ifdef DEBUG
printf("\t%d: %s\n", pl->ngt, ops[pl->cur_token]) ;
#endif
			term(pl, level+1) ;
			add_code(pl, JNZ, level) ;
			break ;
		case ERR : return ;
		default :
			chop_last(pl, JNZ, level) ;
			if(full_expr) {
				add_code(pl, LBL, level) ;
			}
			return ;
		}
	}
}

static void term(plvm pl, int level)
/* compiles a term. A term is a primitive, eventually ANDed
with other primitives */
{
	int full_term = 0 ;

	if(pl->cur_token == ERR) return ;
	prim(pl, level+1) ;
	add_code(pl, JZ, level) ;
	for(;;) {
		switch(pl->cur_token) {
		case AND :
			full_term = 1 ;
			get_token(pl) ;
#ifdef DEBUG
printf("\t%d: %s\n", pl->ngt, ops[pl->cur_token]) ;
#endif
			prim(pl, level+1) ;
			add_code(pl, JZ, level) ;
			break ;
		case ERR : return ;
		default :
			chop_last(pl, JZ, level) ;
			if(full_term) {
				add_code(pl, LBL, level) ;
			}
			return ;
		}
	}
}

static void prim(plvm pl, int level)
/* compiles a primitive. A primitive is either a variable name,
the negation of a primitive, or an expression between parenthesis */
{
	switch(pl->cur_token) {
	case ERR : break ;
	case FRAG : 
		add_code(pl, FR, atoi(pl->cur_string)) ;
		get_token(pl) ;
#ifdef DEBUG
printf("\t%d: %s\n", pl->ngt, ops[pl->cur_token]) ;
#endif
		break ;
	case FRNM :
		add_code(pl, FR, flmaxfrag) ;
		get_token(pl) ;
#ifdef DEBUG
printf("\t%d: %s\n", pl->ngt, ops[pl->cur_token]) ;
#endif
		break ;
	case NOT :
		get_token(pl) ;
#ifdef DEBUG
printf("\t%d:%s\n", pl->ngt, ops[pl->cur_token]) ;
#endif
		prim(pl, level+1) ;
		add_code(pl, INV, -1) ;
		break ;
	case LP :
		get_token(pl) ;
#ifdef DEBUG
printf("\t%d: %s\n", pl->ngt, ops[pl->cur_token]) ;
#endif
		expr(pl, level+1) ;
		if(pl->cur_token != RP) pl->cur_token = ERR ;
		get_token(pl) ;
#ifdef DEBUG
printf("\t%d: %s\n", pl->ngt, ops[pl->cur_token]) ;
#endif
		break ;
	case END : break ;
	default : 
		pl->cur_token = ERR ;
		break ;
	}
}

static popcode jump_search(plvm pl)
/* returns the address of first unprocessed JZ or JNZ opcode,
for which the jump field is still ONE. */
{
	popcode p ;

	for(p=pl->root ; 
		!(p==NULL || 
		((p->mnemo==JZ || p->mnemo==JNZ) && (p->jump==ONE))) ;
		p = p->next) ;
#ifdef DEBUG
if(p == NULL) {
printf("No unprocessed jump\n") ;
} else {
printf("Unprocessed jump at %0X\n", p) ;
}
#endif
	return p ;
}


static void jump_proc(plvm pl)
/* sets jump fields to address of opcode to jump to. */
{
	popcode p, q ;
	int a ;

	if(pl->cur_token == ERR) {
		return ;
	}
	while((p = jump_search(pl)) != NULL) {
		a = p->addr ;
/* stops at next label with same number */
		for(q=p ; !(q->mnemo==LBL && q->addr==a) ; q=q->next) ;
/* jumps after the label opcode */
		q=q->next ;
/* jumps after the following label commands */
		for( ; !(q==NULL || q->mnemo!=LBL) ; q=q->next) ;
/* sets the arrival location of the jump */
		p->jump = q ;

#ifdef DEBUG
printf("Jump from %0X to %0X\n", p, q) ;
#endif
	}
}

static popcode frag_search(plvm pl)
/* searches for a variable opcode whose initial value is still unchanged (-1) */
{
	popcode p ;

	for(p = pl->root ;
		!(p == NULL || (p->mnemo==FR && p->value<0)) ;
		p = p->next) ;
#ifdef DEBUG
if(p != NULL) {
	printf("Found an unprocessed variable: ") ;
	print_code(p) ;
} else {
	printf("No unprocessed variable\n") ;
}
#endif
	return p ;
}

static void frag_proc(plvm pl)
/* arranges pvalue in variable opcodes so that all variables with
the same addr field have their pvalue that point 
to the same value field. Also sets fragc et fragv. */
{
	popcode p, q ;
	int a ;
	int *pv ;
	int i ;

	while((p = frag_search(pl)) != NULL) {
		pl->fragc++ ;
		a = p->addr ;
		pv = &(p->value) ;
		p->pvalue = pv ;
/* p->value > 0 means for the moment that the current opcode contains
a fragment name +1. This will be decoded later in this function
to set fragv. */
		p->value = 1+a ;
		for(q=p->next ; q!=NULL ; q=q->next) {
			if(q->mnemo==FR && q->addr==a) {
				q->pvalue = pv ;
				q->value = 0 ;
#ifdef DEBUG
printf("Identical variable found at: ") ;
print_code(q) ;
#endif
			}
		}
#ifdef DEBUG
printf("Variable processed: ") ;
print_code(p) ;
#endif
	}
/* creates the table of variable names */
	pl->fragv = (int *)malloc(pl->fragc * sizeof(int)) ;
	if(pl->fragv == NULL) {
		pl->cur_token = ERR ;
		return ;
	}
/* scans through opcodes to find unique variable names */
	for(p=pl->root,i=0 ; p!=NULL ; p=p->next) {
		a = p->value ;
		if(a > 0) {
/* variable name is v-1, see above */
			(pl->fragv)[i] = a - 1 ;
			i++ ;
		}
	}
}

static int frag_value(popcode p, lvm_eval f)
/* evaluates variable value for popcode p, according to function f */
{
	int *pv ;
	int v ;

	pv = p->pvalue ;
	v = *pv ;
#ifdef DEBUG
printf("Evaluation of: ") ; print_code(p) ;
#endif
	if(v < 0) {
/* no value cached, f must be called now */
		v = (*f)(p->addr) ;
/* any non 0 value is changed to 1 */
		if(v) v = 1 ;
/* v is cached */
		*pv = v ;
#ifdef DEBUG
printf("New value is: %d\n", v) ;
#endif
	} 
#ifdef DEBUG
	else {
printf("Old value was: %d\n", v) ;
	}
#endif
	return v ;
}
