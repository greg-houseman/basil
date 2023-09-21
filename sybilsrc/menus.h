/*--------------------------------------------------------------------
 *    Basil / Sybil:   menus.h  1.1  1 October 1998
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/

#ifndef _menus_h
#define _menus_h
typedef struct _menuitemdata{
    char        *label;         /* the label for the item */
    WidgetClass *class;         /* pushbutton, label, separator... */
    int		bg;		/* button background colour (-1 dflt) */
    int		spacing;	/* button separation (0 dflt) */
    int     fontindex;  /* index into font list */
    void       (*callback)();   /* routine to call; NULL if none */
    XtPointer    callback_data; /* client_data for callback() */
    struct _menuitemdata *submenu; /* pullright menu items, if not NULL */
    int     submenu_data;
} MenuItemData;
#endif
 
