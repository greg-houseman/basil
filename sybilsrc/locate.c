/*--------------------------------------------------------------------
 *    Basil / Sybil:   locate.c  1.1  1 October 1998
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <Xm/MessageB.h>
#include <Xm/ArrowBG.h>
#include <Xm/PushB.h>
#include <Xm/CascadeB.h>
#include <Xm/Form.h>
#include <Xm/RowColumn.h>
#include <Xm/LabelG.h>
#include "cmndefs.h"
#include "types.h"

#define INCR_I    1
#define DECR_I    2
#define INCR_J    3
#define DECR_J    4

Widget label_i, label_j;

void RowcolDone();
extern void UpdateLocate();
extern void ChangeCell();

void LocateDialog( parent,pw )
Widget parent,pw;
{
    static Widget rowcol_dialog;
    Widget arrow_up_i, arrow_dn_i, arrow_up_j, arrow_dn_j, rowcol;
    char buf[MAX_NUM_LEN+1];
    int i=0;
    Dimension wdth;
    Arg args[2];
    void change_value();

    if (!rowcol_dialog) {
        XtSetArg(args[i], XmNtitle,"Set row and column");i++;
        rowcol_dialog = XmCreateMessageDialog(parent,"rowcol_dialog", args,i);
        XtAddCallback(rowcol_dialog, XmNokCallback, RowcolDone, pw);

        /* No help is available... */
        XtUnmanageChild(XmMessageBoxGetChild(rowcol_dialog,
                XmDIALOG_HELP_BUTTON));
        XtVaGetValues(XmMessageBoxGetChild(rowcol_dialog,
                                        XmDIALOG_OK_BUTTON),
                                         XmNwidth,&wdth,NULL);

        rowcol = XtVaCreateManagedWidget("rowcol",
            xmFormWidgetClass, rowcol_dialog,
            XmNwidth, 4*wdth,
            NULL);

        arrow_up_i = XtVaCreateManagedWidget("arrow_up_i",
            xmArrowButtonGadgetClass, rowcol,
            XmNarrowDirection,   XmARROW_UP,
            XmNtopAttachment,    XmATTACH_FORM,
            XmNleftAttachment,    XmATTACH_FORM,
            NULL);
        XtAddCallback(arrow_up_i, XmNarmCallback, change_value, (XtPointer)INCR_I);

        arrow_dn_i = XtVaCreateManagedWidget("arrow_dn_i",
            xmArrowButtonGadgetClass, rowcol,
            XmNarrowDirection,   XmARROW_DOWN,
            XmNtopAttachment,    XmATTACH_FORM,
            XmNleftAttachment,    XmATTACH_WIDGET,
            XmNleftWidget, arrow_up_i,
            NULL);
        XtAddCallback(arrow_dn_i, XmNarmCallback, change_value, (XtPointer)DECR_I);

        sprintf(buf,"%-4d",Plot_info.curr_cell/(Plot_info.max_col+1));
        label_i = XtVaCreateManagedWidget("label_i",
            xmLabelGadgetClass, rowcol,
            XmNlabelString, XmStringCreateSimple(buf),
            XmNtopAttachment,    XmATTACH_FORM,
            XmNleftAttachment,    XmATTACH_WIDGET,
            XmNleftWidget, arrow_dn_i,
            NULL);

        arrow_dn_j = XtVaCreateManagedWidget("arrow_dn_j",
            xmArrowButtonGadgetClass, rowcol,
            XmNarrowDirection,   XmARROW_DOWN,
            XmNtopAttachment,    XmATTACH_FORM,
            XmNrightAttachment,    XmATTACH_FORM,
            NULL);
        XtAddCallback(arrow_dn_j, XmNarmCallback, change_value, (XtPointer)DECR_J);

        arrow_up_j = XtVaCreateManagedWidget("arrow_up_j",
            xmArrowButtonGadgetClass, rowcol,
            XmNarrowDirection,   XmARROW_UP,
            XmNtopAttachment,    XmATTACH_FORM,
            XmNrightAttachment,    XmATTACH_WIDGET,
            XmNrightWidget, arrow_dn_j,
            NULL);
        XtAddCallback(arrow_up_j, XmNarmCallback, change_value, (XtPointer)INCR_J);

        snprintf(buf,5,"%4u",Plot_info.curr_cell%(Plot_info.max_col+1));
        label_j = XtVaCreateManagedWidget("label_j",
            xmLabelGadgetClass, rowcol,
            XmNlabelString, XmStringCreateSimple(buf),
            XmNtopAttachment,    XmATTACH_FORM,
            XmNrightAttachment,    XmATTACH_WIDGET,
            XmNrightWidget, arrow_up_j,
            NULL);
        XtManageChild (rowcol);
    }
    else {
        snprintf(buf,5,"%-4d",Plot_info.curr_cell/(Plot_info.max_col+1));
        XtVaSetValues(label_i,
            XtVaTypedArg, XmNlabelString, XmRString, buf, strlen(buf),
            NULL);
        sprintf(buf,"%4d",Plot_info.curr_cell%(Plot_info.max_col+1));
        XtVaSetValues(label_j,
            XtVaTypedArg, XmNlabelString, XmRString, buf, strlen(buf),
            NULL);
    }
    XtManageChild (rowcol_dialog);
    XtPopup(XtParent(rowcol_dialog), XtGrabNone);
}

/*
 * rowcol ok button callback function
 */
void RowcolDone(w, client_data, call_data)
Widget w;
XtPointer client_data, call_data;
{
    char *str;
    int row,col,prev;
    XmString xmstr;

    XtVaGetValues(label_i,XmNlabelString,&xmstr,NULL);
    XmStringGetLtoR(xmstr,"",&str);
    row = atoi( str );
    XtFree(str);
    XmStringFree(xmstr);
    XtVaGetValues(label_j,XmNlabelString,&xmstr,NULL);
    XmStringGetLtoR(xmstr,"",&str);
    col = atoi( str );
    XtFree(str);
    XmStringFree(xmstr);
    prev = Plot_info.curr_cell;
    Plot_info.curr_cell = row*(Plot_info.max_col+1) + col;
    if (prev!=Plot_info.curr_cell)
        ChangeCell( prev,Plot_info.curr_cell );
}

/*
 * change_value is called each time an arrow button is selected.
 * If the value has reached its maximum or minimum, just return.
 */
void change_value(w,choice,cbs)
Widget w;
XtPointer choice;
XmArrowButtonCallbackStruct *cbs;
{
    char *str,buf[MAX_NUM_LEN+1],format[8];
    int value,incr,max;
    XmString xmstr;
    Widget dum;

    switch((long)choice) {
    case INCR_I:incr = 1;
                max = Plot_info.max_row;
                dum = label_i;
                strcpy(format,"%4d");
                break;
    case DECR_I:incr = -1;
                max = Plot_info.max_row;
                dum = label_i;
                strcpy(format,"%4d");
                break;
    case INCR_J: incr = 1;
                max = Plot_info.max_col;
                dum = label_j;
                strcpy(format,"%-4d");
                break;
    case DECR_J: incr = -1;
                max = Plot_info.max_col;
                dum = label_j;
                strcpy(format,"%-4d");
                break;
    }
    XtVaGetValues(dum,XmNlabelString,&xmstr,NULL);
    XmStringGetLtoR(xmstr,"",&str);
    value = atoi( str );
    XtFree(str);
    XmStringFree(xmstr);
    if (value + incr > max|| value + incr < 0)
        return;
    value += incr;
    sprintf(buf, format, value);
    XtVaSetValues(dum,
        XtVaTypedArg, XmNlabelString, XmRString, buf, strlen(buf),
        NULL);
}
