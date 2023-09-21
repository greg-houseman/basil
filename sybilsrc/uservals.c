/*--------------------------------------------------------------------
 *    Basil / Sybil:   uservals.c  1.1  1 October 1998
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <Xm/DialogS.h>
#include <Xm/PushB.h>
#include <Xm/PanedW.h>
#include <Xm/Form.h>
#include <Xm/RowColumn.h>
#include <Xm/LabelG.h>
#include <Xm/ToggleBG.h>
#include <Xm/SeparatoG.h>
#include <Xm/TextF.h>
#include "types.h"
#include "errnum.h"
#include "error.h"
#include "xerror.h"
#include "plotdefs.h"
#include "globals.h"
#include "routines.h"

#define COLS  8

typedef struct {
    float *numericaldata;
    char *stringdata;
    int *flag;
    int max_items;
    Widget p_u_shell;
    unsigned char validate;
} DataItem;

typedef struct {
    char *label;
    void (*callback)();
    DataItem *data;
} ActionAreaItem;

static Widget Param[6],cntrvals_dialog;

void checkToggleState();
void checkprflToggleState();
void CntrvalsDialog(Widget, float *, int *);
void CreateActionArea();
void DrawPreviewLine(Widget, int, int, int, int);
void previewline();
void PrflvalsDialog();
void PrflptsDialog();
void prfl_toggled();
void toggled(),set_x(),set_y();
void Translate_u_to_win(float, float, int*, int*);
void UserDone(),LblDone();
void UserCancel();
void XmTextSetString( Widget widget, char * value);
/* Pass in labels and corresponding array of values which will be 
   reset if changed ?? */

void CntrvalsDialog( parent,vals,user_answr )
Widget parent;
float *vals;
int *user_answr;
{
    static Widget pane;
    Widget rowcol1,w,toggle_box;
    char buf[30];
    int i=0;
    Dimension wdth;
    static DataItem user_data;
    static ActionAreaItem button_items[] = {
        { "Ok",     UserDone, NULL          },
        { "Cancel", UserCancel, NULL          },
    };

    /*if (!cntrvals_dialog) {*/
        cntrvals_dialog = XtVaCreatePopupShell("cntrvals_dialog",
            xmDialogShellWidgetClass, parent,
            NULL);

        button_items[0].data = &user_data;
        button_items[1].data = &user_data;
        user_data.numericaldata = vals;

        pane = XtVaCreateWidget("pane",
            xmPanedWindowWidgetClass, cntrvals_dialog,
            XmNsashWidth,  1,
            XmNsashHeight, 1,
            NULL);

        rowcol1 = XtVaCreateWidget("rowcol1",
            xmRowColumnWidgetClass, pane,
            XmNpacking,            XmPACK_COLUMN,
            XmNnumColumns,        COLS,
            XmNorientation,        XmHORIZONTAL,
            XmNisAligned,        True,
            XmNentryAlignment,    XmALIGNMENT_END,
            NULL);

        /* create toggle box for setting max,min */
        toggle_box=XtVaCreateManagedWidget("toggle",
            xmRowColumnWidgetClass, pane, NULL );

        w=XtVaCreateManagedWidget("Use data max, min",
            xmToggleButtonGadgetClass, toggle_box,
            /*XmNuserData, (XtPointer)&user_data,*/
            NULL );
        XtAddCallback( w,XmNvalueChangedCallback,toggled,
                                         (XtPointer)&user_data );
        XmToggleButtonGadgetSetState(w,0,0);

        i=0;
        /* create TextField Widgets with Labels */
        sprintf(buf,"%s%-11.4e%s","Min [",vals[CNTR_MIN],"]:");
        XtVaCreateManagedWidget(buf,
            xmLabelGadgetClass, rowcol1, NULL );

        Param[i] = XtVaCreateManagedWidget("min_val",
            xmTextFieldWidgetClass, rowcol1, 
            XmNuserData, (XtPointer)w,
            XmNcolumns, 13, NULL );
        XtAddCallback(Param[i],XmNvalueChangedCallback,checkToggleState,
                                   (XtPointer)&user_data);
        i++;

        sprintf(buf,"%s%-11.4e%s","Max [",vals[CNTR_MAX],"]:");
        XtVaCreateManagedWidget(buf,
            xmLabelGadgetClass, rowcol1, NULL );

        Param[i] = XtVaCreateManagedWidget("max_val",
            xmTextFieldWidgetClass, rowcol1, 
            XmNuserData, (XtPointer)w,
            XmNcolumns, 13, NULL );
        XtAddCallback(Param[i],XmNvalueChangedCallback,checkToggleState,
                                   (XtPointer)&user_data);
        i++;

        XtVaCreateManagedWidget("Scale:",
            xmLabelGadgetClass, rowcol1, NULL );

        Param[i] = XtVaCreateManagedWidget("scale_val",
            xmTextFieldWidgetClass, rowcol1, 
            XmNcolumns, 13, NULL );
        i++;

        /* create TextField Widgets with Labels */
        sprintf(buf,"%s%-11.4e%s","level [",vals[CNTR_LVL],"]:");
        XtVaCreateManagedWidget(buf,
            xmLabelGadgetClass, rowcol1, NULL );

        Param[i] = XtVaCreateManagedWidget("level",
            xmTextFieldWidgetClass, rowcol1, 
            XmNuserData, (XtPointer)w,
            XmNcolumns, 13, NULL );
        XtAddCallback(Param[i],XmNvalueChangedCallback,checkToggleState,
                                   (XtPointer)&user_data);
        i++;

        sprintf(buf,"%s%-11.4e%s","step [",vals[CNTR_STEP],"]:");
        XtVaCreateManagedWidget(buf,
            xmLabelGadgetClass, rowcol1, NULL );

        Param[i] = XtVaCreateManagedWidget("step",
            xmTextFieldWidgetClass, rowcol1, 
            XmNuserData, (XtPointer)w,
            XmNcolumns, 13, NULL );
        XtAddCallback(Param[i],XmNvalueChangedCallback,checkToggleState,
                                   (XtPointer)&user_data);
        i++;

        XtVaCreateManagedWidget("Max contours:",
            xmLabelGadgetClass, rowcol1, NULL );

        Param[i] = XtVaCreateManagedWidget("max_cntrs",
            xmTextFieldWidgetClass, rowcol1,
            XmNcolumns, 13, NULL );
        i++;

        XtManageChild (toggle_box);
        XtManageChild (rowcol1);

        XtVaGetValues( rowcol1,XmNwidth,&wdth,NULL );

    /* Create the action area */

        CreateActionArea(pane,button_items,XtNumber(button_items),wdth);
        strcpy(buf, "Values used for contouring");
        XtVaSetValues(cntrvals_dialog, 
               XmNtitle, buf,
               NULL);
    /*}*/
    XtManageChild (pane);

    /* display actual data max and min, not user input */
    i=0;
    sprintf(buf,"%11.4e",vals[CNTR_MIN_U]);
    XtVaSetValues(Param[i],
           XmNvalue, buf,
           NULL);

    i++;
    sprintf(buf,"%11.4e",vals[CNTR_MAX_U]);
    XtVaSetValues(Param[i],
           XmNvalue, buf,
           NULL);

    i++;
    sprintf(buf,"%8.4f",vals[CNTR_SCL]);
    XtVaSetValues(Param[i],
           XmNvalue, buf,
           NULL);
    i++;

    sprintf(buf,"%11.4e",vals[CNTR_LVL_U]);
    XtVaSetValues(Param[i],
           XmNvalue, buf,
           NULL);

    i++;
    if (vals[CNTR_STEP_U] > 1000 || vals[CNTR_STEP_U] < 0.001) 
        sprintf(buf,"%11.4e",vals[CNTR_STEP_U]);
    else sprintf(buf,"%8.4f",vals[CNTR_STEP_U]);
    XtVaSetValues(Param[i],
           XmNvalue, buf,
           NULL);

    i++;
    sprintf(buf,"%d",(int)vals[MAX_CNTRS]);
    XtVaSetValues(Param[i],
           XmNvalue, buf,
           NULL);
    i++;

    user_data.flag = user_answr;
    user_data.max_items = i;
    user_data.p_u_shell = cntrvals_dialog;
    user_data.validate = 1;
    XtPopup(cntrvals_dialog, XtGrabNone);
}

void PrflptsDialog( parent,vals,user_answr,draw_line )
Widget parent;
float *vals;
int *user_answr;
unsigned char draw_line;
{
    Widget rowcol1,pane,prfl_dialog;
    char buf[30];
    int i,num_params,num_buttons;
    Dimension wdth;
    ActionAreaItem *buttons;
    static DataItem user_data;
    static ActionAreaItem line_button_items[] = {
        { "Ok",     UserDone, &user_data          },
        { "Preview",     previewline, NULL          },
        { "Cancel", UserCancel, &user_data          },
    };
    static ActionAreaItem button_items[] = {
        { "Ok",     UserDone, &user_data          },
        { "Cancel", UserCancel, &user_data          },
    };

        prfl_dialog = XtVaCreatePopupShell("prfl_dialog",
            xmDialogShellWidgetClass, parent,
            NULL);

        if (draw_line) {
            buttons = line_button_items;
            num_buttons = XtNumber(line_button_items);
        }
        else{
            buttons = button_items;
            num_buttons = XtNumber(button_items);
        }

        pane = XtVaCreateWidget("pane",
            xmPanedWindowWidgetClass, prfl_dialog,
            XmNsashWidth,  1,
            XmNsashHeight, 1,
            NULL);

        rowcol1 = XtVaCreateWidget("rowcol1",
            xmRowColumnWidgetClass, pane,
            XmNpacking,            XmPACK_COLUMN,
            XmNnumColumns,        4,
            XmNorientation,        XmHORIZONTAL,
            XmNisAligned,        True,
            XmNentryAlignment,    XmALIGNMENT_END,
            XmNmarginWidth,        50,
            NULL);


        i=0;
        /* create TextField Widgets with Labels */
        XtVaCreateManagedWidget("       X1:",
            xmLabelGadgetClass, rowcol1, NULL );

        Param[i] = XtVaCreateManagedWidget("x1",
            xmTextFieldWidgetClass, rowcol1, 
            XmNcolumns, 8, NULL );
        i++;

        XtVaCreateManagedWidget("       Y1:",
            xmLabelGadgetClass, rowcol1, NULL );

        Param[i] = XtVaCreateManagedWidget("y1",
            xmTextFieldWidgetClass, rowcol1, 
            XmNcolumns, 8, NULL );
        i++;

        XtVaCreateManagedWidget("       X2:",
            xmLabelGadgetClass, rowcol1, NULL );

        Param[i] = XtVaCreateManagedWidget("x2",
            xmTextFieldWidgetClass, rowcol1, 
            XmNcolumns, 8, NULL );
        i++;

        XtVaCreateManagedWidget("       Y2:",
            xmLabelGadgetClass, rowcol1, NULL );

        Param[i] = XtVaCreateManagedWidget("y2",
            xmTextFieldWidgetClass, rowcol1, 
            XmNcolumns, 8, NULL );
        i++;
        num_params=i;
        XtManageChild (rowcol1);

        XtVaGetValues( rowcol1,XmNwidth,&wdth,NULL );
    /* Create the action area */
           CreateActionArea(pane, buttons, num_buttons, wdth);
        XtManageChild (pane);

    strcpy(buf, "Enter endpoints for profile");
    XtVaSetValues(prfl_dialog, 
           XmNtitle, buf,
           NULL);

    for (i=0;i<num_params;i++) {
         sprintf(buf,"%8.4f",vals[i]);
         XtVaSetValues(Param[i],
                XmNvalue, buf,
                NULL);
    }

    user_data.numericaldata = vals;
    user_data.flag = user_answr;
    user_data.max_items = num_params;
    user_data.p_u_shell = prfl_dialog;
    user_data.validate = 0;
    XtPopup(prfl_dialog, XtGrabNone);
}

void ZoomValsDialog( parent,vals,user_answr )
Widget parent;
float *vals;
int *user_answr;
{
    Widget rowcol1,pane,zoom_dialog;
    char buf[30];
    int i,num_params,num_buttons;
    Dimension wdth;
    ActionAreaItem *buttons;
    static DataItem user_data;
    static ActionAreaItem button_items[] = {
        { "Ok",     UserDone, &user_data          },
        { "Cancel", UserCancel, &user_data          },
    };

        zoom_dialog = XtVaCreatePopupShell("zoom_dialog",
            xmDialogShellWidgetClass, parent,
            NULL);

        buttons = button_items;
        num_buttons = XtNumber(button_items);
        pane = XtVaCreateWidget("pane",
            xmPanedWindowWidgetClass, zoom_dialog,
            XmNsashWidth,  1,
            XmNsashHeight, 1,
            NULL);

        rowcol1 = XtVaCreateWidget("rowcol1",
            xmRowColumnWidgetClass, pane,
            XmNpacking,            XmPACK_COLUMN,
            XmNnumColumns,        4,
            XmNorientation,        XmHORIZONTAL,
            XmNisAligned,        True,
            XmNentryAlignment,    XmALIGNMENT_END,
            XmNmarginWidth,        50,
            NULL);

        i=0;
        /* create TextField Widgets with Labels */
        XtVaCreateManagedWidget("zoom factor:",
            xmLabelGadgetClass, rowcol1, NULL );

        Param[i] = XtVaCreateManagedWidget("x1",
            xmTextFieldWidgetClass, rowcol1, 
            XmNcolumns, 8, NULL );
        i++;

        XtVaCreateManagedWidget(" centre X:",
            xmLabelGadgetClass, rowcol1, NULL );

        Param[i] = XtVaCreateManagedWidget("y1",
            xmTextFieldWidgetClass, rowcol1, 
            XmNcolumns, 8, NULL );
        i++;

        XtVaCreateManagedWidget(" centre Y:",
            xmLabelGadgetClass, rowcol1, NULL );

        Param[i] = XtVaCreateManagedWidget("x2",
            xmTextFieldWidgetClass, rowcol1, 
            XmNcolumns, 8, NULL );
        i++;

        num_params=i;
        XtManageChild (rowcol1);

        XtVaGetValues( rowcol1,XmNwidth,&wdth,NULL );
    /* Create the action area */
        CreateActionArea(pane, buttons, num_buttons, wdth);
        XtManageChild (pane);

       strcpy(buf, "zoom factor / centre");
       XtVaSetValues(zoom_dialog, 
           XmNtitle, buf, NULL);

    for (i=0;i<num_params;i++) {
         sprintf(buf,"%8.4f",vals[i]);
         XtVaSetValues(Param[i],
                XmNvalue, buf, NULL);
    }

    user_data.numericaldata = vals;
    user_data.flag = user_answr;
    user_data.max_items = num_params;
    user_data.p_u_shell = zoom_dialog;
    user_data.validate = 0;
    XtPopup(zoom_dialog, XtGrabNone);
}

void PrflvalsDialog( parent,vals,user_answr )
Widget parent;
float *vals;
int *user_answr;
{
    Widget rowcol1,w,pane,prfl_dialog,toggle_box;
    char buf[30];
    int i,num_params;
    Dimension wdth;
    static DataItem user_data;
    static ActionAreaItem button_items[] = {
        { "Ok",     UserDone, NULL          },
        { "Cancel", UserCancel, NULL          },
    };

        prfl_dialog = XtVaCreatePopupShell("prfl_dialog",
            xmDialogShellWidgetClass, parent,
            NULL);

        button_items[0].data = &user_data;
        button_items[1].data = &user_data;
        user_data.numericaldata = vals;

        pane = XtVaCreateWidget("pane",
            xmPanedWindowWidgetClass, prfl_dialog,
            XmNsashWidth,  1,
            XmNsashHeight, 1,
            NULL);

        rowcol1 = XtVaCreateWidget("rowcol1",
            xmRowColumnWidgetClass, pane,
            XmNpacking,            XmPACK_COLUMN,
            XmNnumColumns,        4,
            XmNorientation,        XmHORIZONTAL,
            XmNisAligned,        True,
            XmNentryAlignment,    XmALIGNMENT_END,
            NULL);

        /* create toggle box for setting max,min */
        toggle_box=XtVaCreateManagedWidget("toggle",
						            xmRowColumnWidgetClass, pane, NULL );

        w=XtVaCreateManagedWidget("Use data max, min",
			            xmToggleButtonGadgetClass, toggle_box,
			           NULL );
        XtAddCallback( w,XmNvalueChangedCallback,prfl_toggled,
                                   (XtPointer)&user_data );
        XmToggleButtonGadgetSetState(w,0,0);

        i=0;
        /* create TextField Widgets with Labels */
        sprintf(buf,"%s%11.4e%s","Min [",vals[PRFL_MIN],"]:");
        XtVaCreateManagedWidget(buf,
            xmLabelGadgetClass, rowcol1, NULL );

        Param[i] = XtVaCreateManagedWidget("min",
            xmTextFieldWidgetClass, rowcol1, 
			XmNuserData, (XtPointer)w,
            XmNcolumns, 11, NULL );
        XtAddCallback(Param[i],XmNvalueChangedCallback,checkprflToggleState,
                                   (XtPointer)&user_data);
        i++;

        sprintf(buf,"%s%11.4e%s","Max [",vals[PRFL_MAX],"]:");
        XtVaCreateManagedWidget(buf,
            xmLabelGadgetClass, rowcol1, NULL );

        Param[i] = XtVaCreateManagedWidget("max",
            xmTextFieldWidgetClass, rowcol1, 
			XmNuserData, (XtPointer)w,
            XmNcolumns, 11, NULL );
        XtAddCallback(Param[i],XmNvalueChangedCallback,checkprflToggleState,
                                   (XtPointer)&user_data);
        i++;

        XtVaCreateManagedWidget("SCALE:",
            xmLabelGadgetClass, rowcol1, NULL );

        Param[i] = XtVaCreateManagedWidget("scale",
            xmTextFieldWidgetClass, rowcol1, 
            XmNcolumns, 11, NULL );
        i++;

        XtVaCreateManagedWidget("SHIFT:",
            xmLabelGadgetClass, rowcol1, NULL );

        Param[i] = XtVaCreateManagedWidget("shift",
            xmTextFieldWidgetClass, rowcol1, 
            XmNcolumns, 11, NULL );
        i++;

    num_params = i;
        XtManageChild (rowcol1);

        XtVaGetValues( rowcol1,XmNwidth,&wdth,NULL );

    /* Create the action area */
        CreateActionArea(pane, button_items, XtNumber(button_items), wdth);
        XtManageChild (pane);

    strcpy(buf, "Enter values for profile");
    XtVaSetValues(prfl_dialog, 
           XmNtitle, buf,
           NULL);

    i=0;
    sprintf(buf,"%11.4e",vals[PRFL_MIN_U]);
    XtVaSetValues(Param[i],
           XmNvalue, buf,
           NULL);
    i++;
    sprintf(buf,"%11.4e",vals[PRFL_MAX_U]);
    XtVaSetValues(Param[i],
           XmNvalue, buf,
           NULL);
    i++;
    sprintf(buf,"%6.4f",vals[i]);
    XtVaSetValues(Param[i],
           XmNvalue, buf,
           NULL);
    i++;
    sprintf(buf,"%6.4f",vals[i]);
    XtVaSetValues(Param[i],
           XmNvalue, buf,
           NULL);
    i++;

    user_data.flag = user_answr;
    user_data.max_items = num_params;
    user_data.p_u_shell = prfl_dialog;
    user_data.validate = 1;
    XtPopup(prfl_dialog, XtGrabNone);
}

Widget label_min, label_max;

void Prfl2D_Dialog( parent,vals,plot,user_answr )
Widget parent;
float *vals;
int *plot;
int *user_answr;
{
    Widget rowcol1,pane,radio_box,radio_x,radio_y;
    char buf[30];
    int i=0;
    Dimension wdth;
    static DataItem user_data;
    static ActionAreaItem button_items[] = {
        { "Ok",     UserDone, NULL          },
        { "Cancel", UserCancel, NULL          },
    };

/*
    if (!cntrvals_dialog) {
*/
        cntrvals_dialog = XtVaCreatePopupShell("cntrvals_dialog",
            xmDialogShellWidgetClass, parent,
            NULL);

        button_items[0].data = &user_data;
        button_items[1].data = &user_data;

        pane = XtVaCreateWidget("pane",
            xmPanedWindowWidgetClass, cntrvals_dialog,
            XmNsashWidth,  1,
            XmNsashHeight, 1,
            NULL);

        rowcol1 = XtVaCreateWidget("rowcol1",
            xmRowColumnWidgetClass, pane,
            XmNpacking,            XmPACK_COLUMN,
            XmNnumColumns,        2,
            XmNorientation,        XmHORIZONTAL,
            XmNisAligned,        True,
            XmNentryAlignment,    XmALIGNMENT_END,
            NULL);

        i=0;
        /* create TextField Widgets with Labels */

        if (*plot==DIM2_X) sprintf(buf,"%s%11.4e%s","Lower limit [",
                                           /*   vals[PRFL_X1],"]:"); */
                                              Pwindo[XCMIN],"]:");
        else sprintf(buf,"%s%11.4e%s","Lower limit [", 
                                          /*     vals[PRFL_Y1],"]:"); */
                                              Pwindo[YCMIN],"]:");
        label_min = XtVaCreateManagedWidget("label_min",
            xmLabelGadgetClass, rowcol1,
            XtVaTypedArg, XmNlabelString, XmRString, buf, strlen(buf)+1,
            NULL );

        Param[i] = XtVaCreateManagedWidget("min_val",
            xmTextFieldWidgetClass, rowcol1, 
            XmNcolumns, 11, NULL );
        i++;

        if (*plot==DIM2_X) sprintf(buf,"%s%11.4e%s","Upper limit [",
                                           /*   vals[PRFL_X2],"]:"); */
                                              Pwindo[XCMAX],"]:");
        else sprintf(buf,"%s%11.4e%s","Upper limit [",
                                          /*   vals[PRFL_Y2],"]:"); */
                                              Pwindo[YCMAX],"]:");
        label_max = XtVaCreateManagedWidget("label_max",
            xmLabelGadgetClass, rowcol1,
            XtVaTypedArg, XmNlabelString, XmRString, buf, strlen(buf)+1,
            NULL );

        Param[i] = XtVaCreateManagedWidget("max_val",
            xmTextFieldWidgetClass, rowcol1, 
            XmNcolumns, 11, NULL );
        i++;

        radio_box = XmCreateRadioBox( pane, "radio_box",NULL,0 );
        XtVaSetValues(radio_box,
           XmNmarginWidth, 40,
           XmNorientation, XmHORIZONTAL,
           NULL);

        radio_x=XtVaCreateManagedWidget("X direction",
            xmToggleButtonGadgetClass, radio_box, NULL );
        XtAddCallback( radio_x,XmNvalueChangedCallback,set_x,plot );

        radio_y=XtVaCreateManagedWidget("Y direction",
            xmToggleButtonGadgetClass, radio_box,
            NULL );
        XtAddCallback( radio_y,XmNvalueChangedCallback,set_y,plot );

        XtManageChild (rowcol1);
        XtManageChild (radio_box);

        XtVaGetValues( rowcol1,XmNwidth,&wdth,NULL );

    /* Create the action area */
           CreateActionArea(pane,button_items,XtNumber(button_items),wdth);
        XtManageChild (pane);
/*
    }
*/
    /* display actual data max and min, not user input */
    strcpy(buf, "Values used for 2D profiles");
    XtVaSetValues(cntrvals_dialog, 
           XmNtitle, buf,
           NULL);
    i=0;
    sprintf(buf,"%11.4e",vals[LWR_LIM]);
    XtVaSetValues(Param[i],
           XmNvalue, buf,
           NULL);
    i++;
    sprintf(buf,"%11.4e",vals[UPR_LIM]);
    XtVaSetValues(Param[i],
           XmNvalue, buf,
           NULL);
    i++;
 
    XmToggleButtonGadgetSetState(radio_x, (*plot==DIM2_X), 0);
    XtVaSetValues(radio_x,
           XmNuserData, vals,
           NULL);
 
    XmToggleButtonGadgetSetState(radio_y, (*plot==DIM2_Y), 0);
    XtVaSetValues(radio_y,
           XmNuserData, vals,
           NULL);
    
    user_data.numericaldata = &vals[LWR_LIM];
    user_data.flag = user_answr;
    user_data.max_items = i;
    user_data.p_u_shell = cntrvals_dialog;
    user_data.validate = 1;
    XtPopup(cntrvals_dialog, XtGrabNone);
}

/*
 * dialog ok button callback function for contour values
 */
void UserDone(w, client_data, call_data)
Widget w;
XtPointer client_data, call_data;
{
    char *str;
    int i;
    float *vals;
    DataItem *data;

    data = (DataItem *)client_data;
    vals = data->numericaldata;
    for (i=0;i<data->max_items;i++) {
        str = XmTextFieldGetString(Param[i]);
        vals[i] = (float)atof( str );
        XtFree(str);
    }

/*
    XtPopdown(data->p_u_shell);
*/
    /*
     * min > max enables reversal of colour scale
     */
/*  if (data->validate && vals[0] >= vals[1])
        ErrorMsg("Max not less than Min",data->p_u_shell,0);
    else {  */
        XtDestroyWidget(data->p_u_shell);
        *(data->flag) = 1;
 /* }  */
}

/*
 * dialog cancel button callback function
 */
void UserCancel(w, client_data, call_data)
Widget w;
XtPointer client_data, call_data;
{
    DataItem *data;

    data = (DataItem *)client_data;
    *(data->flag)=USER_CANCEL;
/*
    XtPopdown(data->p_u_shell);
*/
    XtDestroyWidget(data->p_u_shell);
}

/*
 * cntr dialog toggle button callback function
 */
void checkToggleState(w, client_data, call_data)
Widget w;
XtPointer client_data, call_data;
{
    Widget togglebutton;
    XtPointer ptr;
    char *str,buf[30];
    int min=0,max=1;
    float val1,val2,*vals;
    float limit1,limit2;
    DataItem *data;

    data = (DataItem *)client_data;
    vals = data->numericaldata;

    XtVaGetValues(w,XmNuserData,&ptr,NULL);
    togglebutton = (Widget)ptr;
    str = XmTextFieldGetString(Param[min]);
    val1 = (float)atof( str );
    str = XmTextFieldGetString(Param[max]);
    val2 = (float)atof( str );
    XtFree(str);
    sprintf(buf,"%11.4e",vals[CNTR_MIN]);
    limit1 = atof(buf);
    sprintf(buf,"%11.4e",vals[CNTR_MAX]);
    limit2 = atof(buf);
    if (val1 != limit1 || val2 != limit2)
        XmToggleButtonGadgetSetState(togglebutton,0,0);
    else
        XmToggleButtonGadgetSetState(togglebutton,1,0);
}

/*
 * cntr dialog toggle button callback function
 */
void toggled(w, client_data, call_data)
Widget w;
XtPointer client_data, call_data;
{
    /*
     * assumes the min is first widget and max is second
     * should match widget resource (eg name)
     */
    char buf[30];
    int min=CNTR_MIN_U,max=CNTR_MAX_U;
    int lvl=CNTR_LVL_U,stp=CNTR_STEP_U;
    float *vals;
    DataItem *udata;

    udata = (DataItem *)client_data;
    vals = udata->numericaldata;

    if (XmToggleButtonGadgetGetState(w)) {
        sprintf(buf,"%11.4e",vals[CNTR_MIN]);
        XmTextSetString(Param[min],buf);
        sprintf(buf,"%11.4e",vals[CNTR_MAX]);
        XmTextSetString(Param[max],buf);
        sprintf(buf,"%11.4e",vals[CNTR_LVL]);
        XmTextSetString(Param[lvl],buf);
        sprintf(buf,"%11.4e",vals[CNTR_STEP]);
        XmTextSetString(Param[stp],buf);
    }
    else {
        sprintf(buf,"%11.4e",vals[CNTR_MIN_U]);
        XmTextSetString(Param[min],buf);
        sprintf(buf,"%11.4e",vals[CNTR_MAX_U]);
        XmTextSetString(Param[max],buf);
        sprintf(buf,"%11.4e",vals[CNTR_LVL_U]);
        XmTextSetString(Param[lvl],buf);
        sprintf(buf,"%11.4e",vals[CNTR_STEP_U]);
        XmTextSetString(Param[stp],buf);
    }
}

/*
 * prfl dialog toggle button callback function
 */
void checkprflToggleState(w, client_data, call_data)
Widget w;
XtPointer client_data, call_data;
{
    Widget togglebutton;
    XtPointer ptr;
    char *str,buf[30];
    int min=0,max=1;
    float val1,val2,*vals;
    float limit1,limit2;
    DataItem *data;

    data = (DataItem *)client_data;
    vals = data->numericaldata;

    XtVaGetValues(w,XmNuserData,&ptr,NULL);
    togglebutton = (Widget)ptr;
    str = XmTextFieldGetString(Param[min]);
    val1 = (float)atof( str );
    str = XmTextFieldGetString(Param[max]);
    val2 = (float)atof( str );
    XtFree(str);
    sprintf(buf,"%11.4e",vals[PRFL_MIN]);
    limit1 = atof(buf);
    sprintf(buf,"%11.4e",vals[PRFL_MAX]);
    limit2 = atof(buf);
    if (val1 != limit1 || val2 != limit2)
        XmToggleButtonGadgetSetState(togglebutton,0,0);
    else
        XmToggleButtonGadgetSetState(togglebutton,1,0);
}

/*
 * prfl dialog toggle button callback function
 */
void prfl_toggled(w, client_data, call_data)
Widget w;
XtPointer client_data, call_data;
{
    /*
     * assumes the min is first widget and max is second
     * should match widget resource (eg name)
     */
    char buf[30];
    int max=1,min=0;
    float *vals;
    DataItem *udata;

    udata = (DataItem *)client_data;
    vals = udata->numericaldata;

    if (XmToggleButtonGadgetGetState(w)) {
        sprintf(buf,"%11.4e",vals[PRFL_MIN]);
        XmTextSetString(Param[min],buf);
        sprintf(buf,"%11.4e",vals[PRFL_MAX]);
        XmTextSetString(Param[max],buf);
    }
    else {
        sprintf(buf,"%11.4e",vals[PRFL_MIN_U]);
        XmTextSetString(Param[min],buf);
        sprintf(buf,"%11.4e",vals[PRFL_MAX_U]);
        XmTextSetString(Param[max],buf);
    }
}

/*
 * dialog x-direction radio button callback function
 */
void set_x(w, client_data, call_data)
Widget w;
XtPointer client_data, call_data;
{
    char buf[30];
    int *val = (int *)client_data;
    float *limits;
    XtPointer ptr;

    XmToggleButtonCallbackStruct *data = 
        (XmToggleButtonCallbackStruct *) call_data;
    XtVaGetValues(w,XmNuserData,&ptr,NULL);
    limits = (float *)ptr;
    if (data->set) {
        *val = DIM2_X;
    /*  sprintf(buf,"%s%11.4e%s","Lower limit [", limits[PRFL_MIN_U],"]:"); */
        sprintf(buf,"%s%11.4e%s","Lower limit [", Pwindo[XCMIN],"]:");
        XtVaSetValues(label_min,
            XtVaTypedArg, XmNlabelString, XmRString, buf, strlen(buf)+1,
            NULL );
    /*  sprintf(buf,"%s%11.4e%s","Upper limit [", limits[PRFL_MAX_U],"]:"); */
        sprintf(buf,"%s%11.4e%s","Upper limit [", Pwindo[XCMAX],"]:");
        XtVaSetValues(label_max,
            XtVaTypedArg, XmNlabelString, XmRString, buf, strlen(buf)+1,
            NULL );
    }
}

/*
 * dialog y-direction radio button callback function
 */
void set_y(w, client_data, call_data)
Widget w;
XtPointer client_data, call_data;
{
    char buf[30];
    int *val = (int *)client_data;
    float *limits;
    XtPointer ptr;

    XmToggleButtonCallbackStruct *data = 
        (XmToggleButtonCallbackStruct *) call_data;
    XtVaGetValues(w,XmNuserData,&ptr,NULL);
    limits = (float *)ptr;
    if (data->set) {
        *val = DIM2_Y;
  /*    sprintf(buf,"%s%11.4e%s","Lower limit [",limits[PRFL_MIN],"]:"); */
        sprintf(buf,"%s%11.4e%s","Lower limit [",Pwindo[YCMIN],"]:");
        XtVaSetValues(label_min,
            XtVaTypedArg, XmNlabelString, XmRString, buf, strlen(buf)+1,
            NULL );
  /*    sprintf(buf,"%s%11.4e%s","Upper limit [", limits[PRFL_MAX],"]:");  */
        sprintf(buf,"%s%11.4e%s","Upper limit [", Pwindo[YCMAX],"]:");
        XtVaSetValues(label_max,
            XtVaTypedArg, XmNlabelString, XmRString, buf, strlen(buf)+1,
            NULL );
    }
}

/*
 * dialog preview button callback function
 */
void previewline(w, client_data, call_data)
Widget w;
XtPointer client_data, call_data;
{
    char *str;
    int i,x1,x2,y1,y2;
    int len=1, err=0;
    float xyvals[4];
    float xref,yref;
/*    XREFM and YREFM usage depends on context */
    xref=Data_vars_fl[XREFM];
    yref=Data_vars_fl[YREFM];
    if (Data_vars_int[IDEFTYP]>110) {
       xref=0.;
       yref=0.;
    }

    sybflush_();
    for (i=0;i<4;i++) {
        str = XmTextFieldGetString(Param[i]);
        xyvals[i] = (float)atof( str );
        XtFree(str);
    }
    if (Data_vars_int[NCOMP]==-1) {
        projectxy_(&xyvals[0],&xyvals[1],&xref,&yref,
                   &len,&Data_vars_int[NCOMP],&err);
        projectxy_(&xyvals[2],&xyvals[3],&xref,&yref,
                   &len,&Data_vars_int[NCOMP],&err);
    }
    Translate_u_to_win(xyvals[0],xyvals[1],&x1,&y1);
    Translate_u_to_win(xyvals[2],xyvals[3],&x2,&y2);
    DrawPreviewLine(w,x1,y1,x2,y2);
}

void CreateActionArea(parent, actions, num_actions, wdth)
Widget parent;
ActionAreaItem *actions;
int num_actions;
Dimension wdth;
{
    Widget action_area, widget;
    int i;

    action_area=XtVaCreateWidget("action_area", xmFormWidgetClass, parent,
        XmNfractionBase, num_actions*2 + 1,
        XmNwidth, wdth,
        NULL);

    for (i = 0; i < num_actions; i++) {
        widget = XtVaCreateManagedWidget(actions[i].label,
            xmPushButtonWidgetClass, action_area,
            XmNleftAttachment,       XmATTACH_POSITION,
            XmNleftPosition,         i*2+1,
            XmNtopAttachment,        XmATTACH_FORM,
            XmNbottomAttachment,     XmATTACH_FORM,
/*
            XmNrightAttachment,         XmATTACH_POSITION,
            XmNrightPosition,        i*2 + 2,
*/
            XmNshowAsDefault,        i == 0,
            XmNdefaultButtonShadowThickness, 1,
            NULL);
        if (actions[i].callback)
            XtAddCallback(widget, XmNactivateCallback,
                actions[i].callback, actions[i].data);
        if (i == 0) {
            /* Set the action_area's default button to the first widget
             * created (or, make the index a parameter to the function
             * or have it be part of the data structure). Also, set the
             * pane window constraint for max and min heights so this
             * particular pane in the PanedWindow is not resizable.
             */
            Dimension height, h;
            XtVaGetValues(action_area, XmNmarginHeight, &h, NULL);
            XtVaGetValues(widget, XmNheight, &height, NULL);
            height += 2 * h;
            XtVaSetValues(action_area,
                XmNdefaultButton, widget,
                XmNpaneMaximum,   height,
                XmNpaneMinimum,   height,
                NULL);
        }
    }

    XtManageChild(action_area);
}

#if XY
Widget Lparam[2],lblvals_dialog;
LblvalsDialog( parent,labels,user_answr )
Widget parent;
char *labels;
int *user_answr;
{
    Widget rowcol1,pane;
    char buf[30];
    int i=0;
    Dimension wdth;
    static DataItem user_data;
    static ActionAreaItem button_items[] = {
        { "Ok",     LblDone, NULL          },
        { "Cancel", UserCancel, NULL          },
    };


    if (!lblvals_dialog) {
        lblvals_dialog = XtVaCreatePopupShell("lblvals_dialog",
            xmDialogShellWidgetClass, parent,
            NULL);

        button_items[0].data = &user_data;
        button_items[1].data = &user_data;

        pane = XtVaCreateWidget("pane",
            xmPanedWindowWidgetClass, lblvals_dialog,
            XmNsashWidth,  1,
            XmNsashHeight, 1,
            NULL);

        rowcol1 = XtVaCreateWidget("rowcol1",
            xmRowColumnWidgetClass, pane,
            XmNpacking,            XmPACK_COLUMN,
            XmNnumColumns,        COLS,
            XmNorientation,        XmHORIZONTAL,
            XmNisAligned,        True,
            XmNentryAlignment,    XmALIGNMENT_END,
            NULL);

        i=0;
        /* create TextField Widgets with Labels */
        XtVaCreateManagedWidget("Left label:",
            xmLabelGadgetClass, rowcol1, NULL );

        param[i] = XtVaCreateManagedWidget("left_label",
            xmTextFieldWidgetClass, rowcol1,
            XmNmaxLength,MAXVARNAME,
            XmNcolumns, 8, NULL );
        i++;

        XtVaCreateManagedWidget("Right label:",
            xmLabelGadgetClass, rowcol1, NULL );

        param[i] = XtVaCreateManagedWidget("right_label",
            xmTextFieldWidgetClass, rowcol1,
            XmNmaxLength,MAXVARNAME,
            XmNcolumns, 5, NULL );
        i++;
        XtManageChild (rowcol1);

        XtVaGetValues( rowcol1,XmNwidth,&wdth );

    /* Create the action area */
           CreateActionArea(pane, button_items, XtNumber(button_items), wdth);
        XtManageChild (pane);
    }
    strcpy(buf, "Labels for default positions");
    XtVaSetValues(lblvals_dialog, 
           XmNtitle, buf,
           NULL);
    i=0;
    XtVaSetValues(param[i],
           XmNvalue, &labels[0][0],
           NULL);
    i++;
    XtVaSetValues(param[i],
           XmNvalue, &labels[1][0],
           NULL);
    i++;
 
    user_data.data = labels;
    user_data.flag = user_answr;
    user_data.p_u_shell = lblvals_dialog;
    user_data.vaidate = 0;
    XtPopup(lblvals_dialog, XtGrabNone);
}
/*
 * dialog ok button callback function for contour values
 */
void LblDone(w, client_data, call_data)
Widget w;
XtPointer client_data, call_data;
{
    char *str;
    int i;
    DataItem *data;
    XmString xmstr;

    data = (DataLblItem *)client_data;
    cntr_vals = data->data;
    i=0;
    str = XmTextFieldGetString(param[i]);
    strncpy(&labels[0][0],str,MAXVARNAME);
    XtFree(str);
    i++;
    str = XmTextFieldGetString(param[i]);
    strncpy(&labels[1][0],str,MAXVARNAME);
    XtFree(str);
    i++;

    *(data->flag) = 1;
    XtPopdown(lblvals_dialog);
}
#endif
