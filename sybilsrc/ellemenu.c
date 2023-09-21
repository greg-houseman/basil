#include <stdio.h>
#include <stdlib.h>

#include <Xm/Xm.h>
#ifdef XSYB
#include <X11/Xlib.h>  /* for type XFontStruct in Plot_info - tidy */
#endif

/*
 * Public include files for widgets
 */
#include <Xm/SelectioB.h>

#include "types.h"
#include "cmndefs.h"
/* only needed for log.h which is only needed for valid_terms */
/*#include "strain.h"*/
/*#include "deform.h"*/

#include "string_utils.h"
/*#include "log.h"*/
#include "elle.h"
#include "errnum.h"
#include "error.h"
#include "xerror.h"

int  file_open(file_data**, char*, int, char*);
extern int GetDirSpec();

typedef struct {
    float x;
    float y;
} Coords;

char Lastellemask[SYB_FILENAME_MAX+1];

static void EFileOpenSelectionOK(), EFileSelectionCancel();
extern void ElleFileChosen();
Widget XmCreateFileSelectionDialog();
Widget XmFileSelectionBoxGetChild();

/*
 * maps a file selection box
 */
void OpenElleDialog(parent,title,file_type)
Widget parent;
char *title;
int file_type;
{
    static Widget fileSelection;

    char *buf;
    XmString tmp1, tmp2;
    Arg args[5];
    int i=0;

    tmp1 = XmStringCreateSimple(title);
    tmp2=NULL;
    XtSetArg(args[i], XmNdialogTitle, tmp1); i++;
    XtSetArg(args[i], XmNuserData, file_type); i++;
    buf = Lastellemask;
    if (!fileSelection) {
      strcpy(Lastellemask,"");
      if (GetDirSpec( buf,file_type )) {
          tmp2 = XmStringCreateLocalized(buf);
          XtSetArg(args[i], XmNdirMask, tmp2); i++;
      }
      fileSelection= (Widget) XmCreateFileSelectionDialog(
            parent,  /* parent widget */
            "fileSelection",  /* widget name */
            args,   /* argument list*/
            i   /* arglist size */
            );
      XtUnmanageChild((Widget) XmFileSelectionBoxGetChild(fileSelection,
                   XmDIALOG_HELP_BUTTON));
      XtAddCallback(fileSelection, XmNcancelCallback,
                    EFileSelectionCancel, NULL);
      XtAddCallback(fileSelection, XmNokCallback, EFileOpenSelectionOK,
                     NULL);
    }
    else {
      if (strlen(buf)==0) GetDirSpec( buf,file_type );
      tmp2 = XmStringCreateLocalized(buf);
      XtSetArg(args[i], XmNdirMask, tmp2); i++;
      XtSetValues(fileSelection,args,i);
    }
    XmStringFree(tmp1); XmStringFree(tmp2);
    XtManageChild(fileSelection);
    XtPopup(XtParent(fileSelection), XtGrabNone);
}


/*
 * EFileSelectionCancel() simply unmaps the file selection box.
 */
static void EFileSelectionCancel(w, client_data, call_data)
Widget w;
XtPointer client_data;
XtPointer call_data;
{
    if (Plot_info.elle_file) Plot_info.elle_file->fp = 0;
    XtUnmanageChild(w);
}       /* EFileSelectionCancel */


static void EFileOpenSelectionOK(w, client_data, cbs)
Widget w;
XtPointer client_data;
XmFileSelectionBoxCallbackStruct *cbs;/*always third param*/
{
    char *filename,*buf;
    int err=0,type=0;
    XtPointer ptr;
    XmString xmstr;

    if (!XmStringGetLtoR(cbs->value, "", &filename))
                return;
    XtVaGetValues(w,XmNuserData,&ptr,NULL);
    type = (long)ptr;
    XtVaGetValues(w,XmNdirMask,&xmstr,NULL);
    XmStringGetLtoR(xmstr,"",&buf);
    XmStringFree(xmstr);
    if (filename[strlen(filename)-1]!='/') {
        err = file_open(&(Plot_info.elle_file),filename,SYB_FILENAME_MAX,"r");
        if (!err) {
            strncpy(Lastellemask,buf,SYB_FILENAME_MAX);
            XtUnmanageChild(w);
            ElleFileChosen();
        }
        else ErrorMsg(filename,w,err);
    }
    else ErrorMsg("The file selected is a directory",w,0);
    XtFree( buf );
    XtFree( filename );
}       /* EFileOpenSelectionOK */

