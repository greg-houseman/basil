/*--------------------------------------------------------------------
 *    Basil / Sybil:   help.c  1.1  1 October 1998
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/

/* use a ScrolledText object to view the
 * contents of a file of help text
 */
#include <X11/Xos.h>
#include <Xm/Text.h>
#include <Xm/TextF.h>
#include <Xm/PanedW.h>
#include <Xm/DialogS.h>
#include <Xm/RowColumn.h>
#include <Xm/Form.h>
#include <Xm/LabelG.h>
#include <Xm/PushB.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include "help.h"
#include "errnum.h"

#define SYB_FILENAME_MAX  1023
#define HELP_FILE "sybilhelp.txt"

void ErrorMsg(char*, Widget, int);

Widget GetTopShell();
void   read_help();
static void DestroyShell();

void HelpView( parent,fg,bg )
Widget parent;
Pixel fg, bg;
{
    char helpfile[SYB_FILENAME_MAX+1],*dir;
    Widget        help_dialog, help_w, button, form, text_w;
    int i=0;
    Arg           args[8];

    if ((dir = getenv("BASILPATH"))==0) {
		ErrorMsg("Please set BASILPATH",parent,OPEN_ERR);
		return;
	}
    /* 
     * ScrolledText/TextField as WorkWindow.
     */
    help_dialog = XtVaCreatePopupShell("Help",
        xmDialogShellWidgetClass, GetTopShell(parent),
        XtNinput, (XtArgVal)True,
        XmNdeleteResponse, XmDESTROY,
        NULL);

   /* MainWindow  - 
     * ScrolledText/TextField as WorkWindow.
     * add search dialog
     */
    help_w = XtVaCreateWidget("help_w",
        xmPanedWindowWidgetClass, help_dialog,
        XmNsashWidth, 1,
        XmNsashHeight, 1,
        XmNforeground, fg,
        XmNbackground, bg,
        /* XmNscrollingPolicy,   XmVARIABLE, */
        NULL);
    form = XtVaCreateManagedWidget("work_area",
        xmFormWidgetClass, help_w,
        XmNforeground, fg,
        XmNbackground, bg,
        NULL);

    /* Create ScrolledText -- this is work area for the MainWindow */
    XtSetArg(args[i], XmNscrollVertical,        True); i++;
    XtSetArg(args[i], XmNscrollHorizontal,      False); i++;
    XtSetArg(args[i], XmNrows,      24); i++;
    XtSetArg(args[i], XmNcolumns,   70); i++;
    XtSetArg(args[i], XmNeditable,  False); i++;
    XtSetArg(args[i], XmNeditMode,  XmMULTI_LINE_EDIT); i++;
    XtSetArg(args[i], XmNforeground,  fg); i++;
    XtSetArg(args[i], XmNbackground,  bg); i++;
    text_w = XmCreateScrolledText(form, "text_w", args, i);

    strcpy(helpfile,dir);
    strcat(helpfile,"/docs/");
    strcat(helpfile,HELP_FILE);
    read_help( text_w,helpfile );

    XtManageChild(text_w);
    XtManageChild(form);
    /* Create another form to act as the action area for the dialog */
    form = XtVaCreateWidget("form2", xmFormWidgetClass, help_w,
        XmNfractionBase,    7,
        XmNforeground, fg,
        XmNbackground, bg,
        NULL);

    /* The Ok button is under the pane's separator and is
     * attached to the left edge of the form.  It spreads from
     * position 1 to 2 along the bottom (the form is split into
     * 7 separate grids by XmNfractionBase ).
     */
    button = XtVaCreateManagedWidget("Ok",
        xmPushButtonWidgetClass, form,
        XmNtopAttachment,        XmATTACH_FORM,
        XmNbottomAttachment,     XmATTACH_FORM,
        XmNleftAttachment,       XmATTACH_POSITION,
        XmNleftPosition,         1,
        XmNrightAttachment,      XmATTACH_POSITION,
        XmNrightPosition,        2,
        XmNshowAsDefault,        True,
        XmNdefaultButtonShadowThickness, 1,
        NULL);
    XtAddCallback(button, XmNactivateCallback, DestroyShell, help_dialog);

    /*
     * This is created with its XmNsensitive resource set to False
     * attach action for searching help to this button
     */
    button = XtVaCreateManagedWidget("Search",
        xmPushButtonWidgetClass, form,
        XmNsensitive,            False,
        XmNtopAttachment,        XmATTACH_FORM,
        XmNbottomAttachment,     XmATTACH_FORM,
        XmNleftAttachment,       XmATTACH_POSITION,
        XmNleftPosition,         3,
        XmNrightAttachment,      XmATTACH_POSITION,
        XmNrightPosition,        4,
        XmNshowAsDefault,        False,
        XmNdefaultButtonShadowThickness, 1,
        NULL);

    /*
     * This is created with its XmNsensitive resource set to False
     * attach action for displaying more help to this button
     */
    button = XtVaCreateManagedWidget("More",
        xmPushButtonWidgetClass, form,
        XmNsensitive,            False,
        XmNtopAttachment,        XmATTACH_FORM,
        XmNbottomAttachment,     XmATTACH_FORM,
        XmNleftAttachment,       XmATTACH_POSITION,
        XmNleftPosition,         5,
        XmNrightAttachment,      XmATTACH_POSITION,
        XmNrightPosition,        6,
        XmNshowAsDefault,        False,
        XmNdefaultButtonShadowThickness, 1,
        NULL);
 
    /* Fix the action area pane to its current height */
    XtManageChild(button);
    XtManageChild(form);
    {
        Dimension h;
        XtVaGetValues(button, XmNheight, &h, NULL);
        XtVaSetValues(form, XmNpaneMaximum, h, XmNpaneMinimum, h, NULL);
   }
    XtManageChild(help_w);
    XtPopup(XtParent(help_dialog),XtGrabNone);
    /*
     *  help window kept appearing behind main window
     */
    XRaiseWindow(XtDisplay(help_dialog),XtWindow(help_dialog));

}

void read_help(text_w, filename)
Widget text_w;
char *filename;
{
    char *text;
    struct stat statb;
    FILE *fp;

    /* make sure the file is a regular text file and open it */
    if (stat(filename, &statb) == -1 ||
            (statb.st_mode & S_IFMT) != S_IFREG ||
            !(fp = fopen(filename, "r"))) {
        if ((statb.st_mode & S_IFMT) == S_IFREG)
            perror(filename); /* send to stderr why we can't read it */
        else
            fprintf(stderr, "%s: not a regular file\n", filename);
        return;
    }

    /* put the contents of the file in the Text widget by allocating
     * enough space for the entire file, reading the file into the
     * allocated space, and using XmTextFieldSetString() to show the file.
     */
    if (!(text = XtMalloc((Cardinal)(statb.st_size+1)))) {
        fprintf(stderr, "Can't alloc enough space for %s", filename);
        fclose(fp);
        return;
    }

    if (!fread(text, sizeof(char), statb.st_size+1, fp))
        fprintf(stderr, "Warning: may not have read entire file!\n");

    text[statb.st_size] = '\0'; /* be sure to NULL-terminate */

    /* insert file contents in Text widget */
    XmTextSetString(text_w, text);

    /* free all allocated space */
    XtFree(text);
    fclose(fp);
}

static void DestroyShell(widget, shell)
Widget widget, shell;
{
    /* fclose(fp); if using more and search buttons */
    XtDestroyWidget(shell);
}

Widget GetTopShell(w)
Widget w;
{
    while ( w && !XtIsWMShell(w) )
    w = XtParent(w);
    return( w );
}
