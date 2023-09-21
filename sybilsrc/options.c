/*--------------------------------------------------------------------
 *    Basil / Sybil:   options.c  1.1  1 October 1998
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
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
#include "globals.h"
#include "plotdefs.h"
#include "errnum.h"
#include "error.h"
#include "strain.h"
#include "deform.h"
#include "string_utils.h"
#include "log.h"
#include "plot.h"

#define COLS  8
#define PARAMS_MAX  16

typedef struct {
    input_options *data;
    int num_ints;
    int num_toggles;
    int num_items;
    Widget p_u_shell;
} OptionDataItem;

typedef struct {
    int id;
    int *initial_addr;
} IntegerOptionItem;

typedef struct {
    int id;
    int *initial_addr;
    int *tmp_addr;
} ToggleItem;

typedef struct {
    int id;
    unsigned char *initial_addr;
    unsigned char *tmp_addr;
} ToggleBoolItem;

typedef struct {
    char *label;
    void (*callback)();
    OptionDataItem *data;
} PlotActionAreaItem;

Widget Opts_params[PARAMS_MAX];
input_options Tmp;

void PlotOptsDialog();
void init_box(float* pwin, input_options* opts);
static void CreateActionArea();
static void PlotDlgDone();
static void UserCancel();
void toggled_int();
void toggled_uchar();
void toggled_lines();
void toggled_shade();
void bar_cb();
void set_tic();
/* void SetParameter(); */

IntegerOptionItem  Option_items1[] = {
    { O_NX3, &Settings.plot_opts.nx3 },
    { O_MP, &Settings.plot_opts.mp },
    { O_MPE, &Settings.plot_opts.mpe },
    { O_NP, &Settings.plot_opts.np },
    };
IntegerOptionItem  Option_items2[] = {
    { O_PRFLPTS, &Settings.plot_opts.profile_pts },
    { O_STIP, &Settings.plot_opts.stipple },
    { O_SOLNROT, &Settings.plot_opts.solution_rot },
    };

ToggleItem Toggle_items1[] = {
    { O_DBLE, &Settings.plot_opts.dble, &Tmp.plot_opts.dble },
    { O_LABEL, &Settings.plot_opts.label, &Tmp.plot_opts.label },
    { O_FLIP, &Settings.plot_opts.flip, &Tmp.plot_opts.flip },
    { O_VERBOSE, &Settings.verbose, &Tmp.verbose },
    };

ToggleBoolItem Toggle_items2[] = {
    { O_CLIPTOCELL, &Settings.clip_to_cell, &Tmp.clip_to_cell },
    { O_MARKCELL, &Settings.mark_cell, &Tmp.mark_cell },
    };

void PlotOptsDialog( w, client_data, call_data )
Widget w;
XtPointer client_data;
XtPointer call_data;
{
    XmString menu_label,labels[6];
    static Widget opts_dialog,bar_menu,tic_menu;
    Widget *parent,rowcol1,rowcol2,pane,toggle_box;
    char buf[30], wdgt_label[MAX_OPTION_NAME+1];
    int i=0, j, state=0;
    int val=0;
    Dimension wdth;
    Position x,y;
    static OptionDataItem user_data;
    static PlotActionAreaItem button_items[] = {
        { "Ok",     PlotDlgDone, NULL          },
        { "Cancel", UserCancel, NULL          },
    };

    Tmp = Settings;
    if (!opts_dialog) {
        parent = (Widget *)client_data;
        opts_dialog = XtVaCreatePopupShell("opts_dialog",
            xmDialogShellWidgetClass, *parent,
            NULL);

        button_items[0].data = &user_data;
        button_items[1].data = &user_data;

        pane = XtVaCreateWidget("pane",
            xmPanedWindowWidgetClass, opts_dialog,
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

        /*
         * create TextField Widgets with Labels
         */
        i=0;
        for (j=0;j<XtNumber(Option_items1);j++) {
            id_match(option_terms,Option_items1[j].id,wdgt_label);
            XtVaCreateManagedWidget(wdgt_label,
                xmLabelGadgetClass, rowcol1, NULL );

            sprintf(buf,"%4d",*(Option_items1[j].initial_addr));
            Opts_params[i] = XtVaCreateManagedWidget("txt_wdgt",
                xmTextFieldWidgetClass, rowcol1, 
                XmNuserData, (XtPointer)&Option_items1[j].id,
                XmNvalue, buf,
                XmNcolumns, 4, NULL );
            i++;
        }

        rowcol2 = XtVaCreateWidget("rowcol2",
            xmRowColumnWidgetClass, pane,
            XmNpacking,            XmPACK_COLUMN,
            XmNnumColumns,        3,
            XmNorientation,        XmHORIZONTAL,
            XmNisAligned,        True,
            XmNentryAlignment,    XmALIGNMENT_END,
            NULL);

        for (j=0;j<XtNumber(Option_items2);j++) {
            id_match(option_terms,Option_items2[j].id,wdgt_label);
            strcat(wdgt_label,":");
            XtVaCreateManagedWidget(wdgt_label,
                xmLabelGadgetClass, rowcol2, NULL );

            sprintf(buf,"%4d",*(Option_items2[j].initial_addr));
            Opts_params[i] = XtVaCreateManagedWidget("txt_wdgt",
                xmTextFieldWidgetClass, rowcol2, 
                XmNuserData, (XtPointer)&Option_items2[j].id,
                XmNvalue, buf,
                XmNcolumns, 4, NULL );
            i++;
        }
/*
 * store the number of integer plot parameters
 */
        user_data.num_ints = i;

        toggle_box=XtVaCreateManagedWidget("toggle",
            xmRowColumnWidgetClass, pane, NULL );

        /*
         * create ToggleButton gadgets with Labels
         */
        for (j=0;j<XtNumber(Toggle_items1);j++) {
            id_match(option_terms,Toggle_items1[j].id,wdgt_label);
            Opts_params[i] = XtVaCreateManagedWidget(wdgt_label,
                xmToggleButtonGadgetClass, toggle_box,
                XmNuserData, (XtPointer)&Toggle_items1[j].id,
                NULL );

            XmToggleButtonGadgetSetState(Opts_params[i],
                           (*(Toggle_items1[j].initial_addr)==1), 0);
            XtAddCallback( Opts_params[i],XmNvalueChangedCallback,
                            toggled_int, Toggle_items1[j].tmp_addr );
            i++;
        }
/*
 * store the number of boolean plot parameters
 */
        user_data.num_toggles = i-user_data.num_ints;

        XtVaCreateManagedWidget("sep1",
            xmSeparatorGadgetClass, toggle_box, NULL );

        XtVaCreateManagedWidget("Contour Options:",
            xmLabelGadgetClass, toggle_box, NULL );

        /*
         * create ToggleButton gadgets for contour plot options
         */
        val = O_CNTRPLT;
        Opts_params[i] = XtVaCreateManagedWidget("Lines",
            xmToggleButtonGadgetClass, toggle_box,
            XmNuserData, (XtPointer)&val,
            NULL );
        XtAddCallback( Opts_params[i],XmNvalueChangedCallback,toggled_lines,
                                         &Tmp.plot_opts.contour_plot );
        XmToggleButtonGadgetSetState(Opts_params[i],
           (Settings.plot_opts.contour_plot&LINES)!=0, 0);
        i++;

        Opts_params[i] = XtVaCreateManagedWidget("Shading",
            xmToggleButtonGadgetClass, toggle_box,
            XmNuserData, (XtPointer)&val,
            NULL );
        XtAddCallback( Opts_params[i],XmNvalueChangedCallback,toggled_shade,
                                         &Tmp.plot_opts.contour_plot );
        XmToggleButtonGadgetSetState(Opts_params[i],
           (Settings.plot_opts.contour_plot&SHADE)!=0, 0);
        i++;

        /*
         * create the option menu for specifying bar orientation
         */
#if XmVersion<1002
        menu_label = XmStringCreateSimple("Bar ");
        for (j=0;j<(XtNumber(colourbar_terms)-1);j++)
          labels[j] = XmStringCreateSimple(colourbar_terms[j].name);
#else
        menu_label = XmStringCreateLocalized("Bar ");
        for (j=0;j<(XtNumber(colourbar_terms)-1);j++)
          labels[j] = XmStringCreateLocalized(colourbar_terms[j].name);
#endif
        if (id_match(colourbar_terms,Settings.plot_opts.colour_bar,buf)) {
          while (state<(XtNumber(colourbar_terms)-1) &&
                      strcmp(buf,colourbar_terms[state].name)) state++;
          if (state==(XtNumber(colourbar_terms)-1)) state = 0;
        }
        bar_menu = XmVaCreateSimpleOptionMenu(toggle_box, "bar_menu",
            menu_label, 'B', state /*initial menu selection*/, bar_cb,
            XmVaPUSHBUTTON, labels[0], 'v', NULL, NULL,
            XmVaPUSHBUTTON, labels[1], 'h', NULL, NULL,
            XmVaPUSHBUTTON, labels[2], 'n', NULL, NULL,
            NULL);
        for (j=0;j<(XtNumber(colourbar_terms)-1);j++) XmStringFree(labels[j]);
        XmStringFree(menu_label);
        XtManageChild(bar_menu);

        /*
         * create the option menu for specifying tic placement
         */
#if XmVersion<1002
        menu_label = XmStringCreateSimple("Tics ");
        for (j=0;j<(XtNumber(ticmark_terms)-1);j++)
          labels[j] = XmStringCreateSimple(ticmark_terms[j].name);
#else
        menu_label = XmStringCreateLocalized("Tics ");
        for (j=0;j<(XtNumber(ticmark_terms)-1);j++)
          labels[j] = XmStringCreateLocalized(ticmark_terms[j].name);
#endif
        if (id_match(ticmark_terms,Settings.plot_opts.ticmark,buf)) {
          while (state<(XtNumber(ticmark_terms)-1) &&
                       strcmp(buf,ticmark_terms[state].name)) state++;
          if (state==(XtNumber(ticmark_terms)-1)) state = 0;
        }
        tic_menu = XmVaCreateSimpleOptionMenu(toggle_box, "tic_menu",
            menu_label, 'T', state /*initial menu selection*/, set_tic,
            XmVaPUSHBUTTON, labels[0], 'c', NULL, NULL,
            XmVaPUSHBUTTON, labels[1], 'i', NULL, NULL,
            XmVaPUSHBUTTON, labels[2], 'e', NULL, NULL,
            XmVaPUSHBUTTON, labels[3], 'n', NULL, NULL,
            NULL);
        for (j=0;j<(XtNumber(ticmark_terms)-1);j++) XmStringFree(labels[j]);
        XmStringFree(menu_label);
        XtManageChild(tic_menu);

        XtVaCreateManagedWidget("sep",
            xmSeparatorGadgetClass, toggle_box, NULL );

        /*
         * create ToggleButton gadgets with Labels
         */
        id_match(option_terms,Toggle_items2[0].id,wdgt_label);
        Opts_params[i] = XtVaCreateManagedWidget(wdgt_label,
                xmToggleButtonGadgetClass, toggle_box,
                XmNuserData, (XtPointer)&Toggle_items2[0].id,
                NULL );

        XmToggleButtonGadgetSetState(Opts_params[i],
                           (*(Toggle_items2[0].initial_addr)==1), 0);
        XtAddCallback( Opts_params[i],XmNvalueChangedCallback,
                            toggled_uchar, Toggle_items2[0].tmp_addr );
        i++;
/*
 * store the number of boolean plot parameters
 */
        user_data.num_toggles++;

        XtVaCreateManagedWidget("sep1",
            xmSeparatorGadgetClass, toggle_box, NULL );

        id_match(option_terms,Toggle_items2[1].id,wdgt_label);
        Opts_params[i] = XtVaCreateManagedWidget(wdgt_label,
                xmToggleButtonGadgetClass, toggle_box,
                XmNuserData, (XtPointer)&Toggle_items2[1].id,
                NULL );

        XmToggleButtonGadgetSetState(Opts_params[i],
                           (*(Toggle_items2[1].initial_addr)==1), 0);
        XtAddCallback( Opts_params[i],XmNvalueChangedCallback,
                            toggled_uchar, Toggle_items2[1].tmp_addr );
        i++;
        user_data.num_toggles++;

        user_data.num_items = user_data.num_ints+user_data.num_toggles;

        XtManageChild (rowcol1);
        XtManageChild (rowcol2);
        XtManageChild (toggle_box);

        XtVaGetValues( rowcol1,XmNwidth,&wdth,NULL );

    /* Create the action area */
        CreateActionArea(pane,button_items,XtNumber(button_items),wdth);

        XtTranslateCoords( w,(Position)0,
                   (Position)0,&x,&y);
        XtVaSetValues( opts_dialog,XmNx,x,
                   XmNy,y,
                   NULL );
        XtManageChild (pane);
        user_data.data = &Tmp;
        user_data.p_u_shell = opts_dialog;
    }
    else {
        if (id_match(colourbar_terms,Settings.plot_opts.colour_bar,buf)) {
#if XmVersion<1002
            labels[0] = XmStringCreateSimple(buf);
#else
            labels[0] = XmStringCreateLocalized(buf);
#endif
            XtVaSetValues( XmOptionButtonGadget(bar_menu),
                           XmNlabelString,labels[0],NULL );
            XmStringFree(labels[0]);
        }
        if (id_match(ticmark_terms,Settings.plot_opts.ticmark,buf)) {
#if XmVersion<1002
            labels[0] = XmStringCreateSimple(buf);
#else
            labels[0] = XmStringCreateLocalized(buf);
#endif
            XtVaSetValues( XmOptionButtonGadget(tic_menu),
                           XmNlabelString,labels[0],NULL );
            XmStringFree(labels[0]);
        }
    }

    XtPopup(opts_dialog, XtGrabNone);
}

/*
 * dialog ok button callback function for plotting options
 * gets the current plotting options and calls NewMesh()
 * if the size of the interpolation mesh has been changed
 */
static void PlotDlgDone(w, client_data, call_data)
Widget w;
XtPointer client_data, call_data;
{
    unsigned char reset = 0;
    char *str;
    int i=0, val, choice;
    OptionDataItem *user_data;
    XtPointer ptr;
    extern Widget plotArea;
    extern void NewMesh(),DrawBorder(),UndrawBorder();

    user_data = (OptionDataItem *)client_data;
    
    /*
     * get values of integer parameters
     * save any changed values for log file
     */
    while (i<PARAMS_MAX && Opts_params[i]!=NULL) {
        XtVaGetValues(Opts_params[i],
               XmNuserData, &ptr,
               NULL);
        choice = *(int *)ptr;
        switch(choice) {
            case O_NX3: str = XmTextFieldGetString(Opts_params[i]);
                        val = atoi(str);
                        XtFree(str);
                        Tmp.plot_opts.nx3 = val;
                        if (val!=Settings.plot_opts.nx3)
                            Save_option((XtPointer)&Tmp.plot_opts.nx3,
                                                    O_NX3,SYB_INT);
                    break;
            case O_MP : str = XmTextFieldGetString(Opts_params[i]);
                        val = atoi(str);
                        XtFree(str);
                        Tmp.plot_opts.mp = val;
                        if (val!=Settings.plot_opts.mp)
                            Save_option((XtPointer)&Tmp.plot_opts.mp,
                                                     O_MP,SYB_INT);
                    break;
            case O_MPE: str = XmTextFieldGetString(Opts_params[i]);
                        val = atoi(str);
                        XtFree(str);
                        Tmp.plot_opts.mpe = val;
                        if (val!=Settings.plot_opts.mpe)
                            Save_option((XtPointer)&Tmp.plot_opts.mpe,
                                                     O_MPE,SYB_INT);
                    break;
            case O_NP : str = XmTextFieldGetString(Opts_params[i]);
                        val = atoi(str);
                        XtFree(str);
                        Tmp.plot_opts.np = val;
                        if (val!=Settings.plot_opts.np)
                            Save_option((XtPointer)&Tmp.plot_opts.np,
                                                     O_NP,SYB_INT);
                    break;
            case O_PRFLPTS : str = XmTextFieldGetString(Opts_params[i]);
                        val = atoi(str);
                        XtFree(str);
                        Tmp.plot_opts.profile_pts = val;
                        if (val!=Settings.plot_opts.profile_pts)
                            Save_option((XtPointer)&Tmp.plot_opts.profile_pts,
                                                     O_PRFLPTS,SYB_INT);
                    break;
            case O_STIP : str = XmTextFieldGetString(Opts_params[i]);
                        val = atoi(str);
                        XtFree(str);
                        Tmp.plot_opts.stipple = val;
                        if (val!=Settings.plot_opts.stipple)
                            Save_option((XtPointer)&Tmp.plot_opts.stipple,
                                                        O_STIP,SYB_INT);
                    break;
            case O_SOLNROT : str = XmTextFieldGetString(Opts_params[i]);
                        val = atoi(str);
                        XtFree(str);
                        Tmp.plot_opts.solution_rot = val;
                        if (val!=Settings.plot_opts.solution_rot) {
                            Save_option((XtPointer)&Tmp.plot_opts.solution_rot,
                                                       O_SOLNROT,SYB_INT);
                            Settings.rescale = 1;
                        }
                    break;
            case O_DBLE : if (Tmp.plot_opts.dble!=Settings.plot_opts.dble)
                            Save_option((XtPointer)&Tmp.plot_opts.dble,
                                                       O_DBLE,SYB_INT);
                    break;
            case O_LABEL : if (Tmp.plot_opts.label!=Settings.plot_opts.label)
                            Save_option((XtPointer)&Tmp.plot_opts.label,
                                                       O_LABEL,SYB_INT);
                    break;
            case O_FLIP : if (Tmp.plot_opts.flip!=Settings.plot_opts.flip)
                            Save_option((XtPointer)&Tmp.plot_opts.flip,
                                                       O_FLIP,SYB_INT);
                    break;
            case O_CLIPTOCELL: if (Tmp.clip_to_cell!=Settings.clip_to_cell) {
                            Save_option((XtPointer)&Tmp.clip_to_cell,
                                                   O_CLIPTOCELL,SYB_BOOL);
                            Settings.clip_to_cell=Tmp.clip_to_cell;
                            Settings.rescale = 0;
                            init_box(Pwindo,&Settings);
                            }
                    break;
            case O_VERBOSE: if (Tmp.verbose!=Settings.verbose)
                            Save_option((XtPointer)&Tmp.verbose,
                                                       O_VERBOSE,SYB_INT);
                    break;
            case O_MARKCELL: if (Tmp.mark_cell!=Settings.mark_cell)
                            Save_option((XtPointer)&Tmp.mark_cell,
                                                       O_MARKCELL,SYB_BOOL);
                    break;
            default: break;
        }
        i++;
    }

    if (Tmp.plot_opts.nx3!=Settings.plot_opts.nx3 ||
        Tmp.plot_opts.dble!=Settings.plot_opts.dble) reset = 1;
    if (Tmp.plot_opts.contour_plot!=Settings.plot_opts.contour_plot) 
        Save_option((XtPointer)&Tmp.plot_opts.contour_plot,O_CNTRPLT,
                                                    SYB_INT);
    if (Tmp.plot_opts.colour_bar!=Settings.plot_opts.colour_bar) 
        Save_option((XtPointer)&Tmp.plot_opts.colour_bar,O_COLBAR,
                                                    SYB_INT);
    if (Tmp.plot_opts.ticmark!=Settings.plot_opts.ticmark) 
        Save_option((XtPointer)&Tmp.plot_opts.ticmark,O_TIC,
                                                    SYB_INT);
    
    Settings = Tmp;

    if (reset) NewMesh();

    if (Settings.mark_cell) DrawCellBorder(plotArea,Plot_info.curr_cell);
    else UndrawCellBorder(plotArea,Plot_info.curr_cell);
    XtPopdown(user_data->p_u_shell);
}

/*
 * dialog cancel button callback function for plotting options
 * sets the current dialog values back to those of Settings
 */
static void UserCancel(w, client_data, call_data)
Widget w;
XtPointer client_data, call_data;
{
    char *str,buf[30];
    int i=0, val, choice;
    OptionDataItem *user_data;
    XtPointer ptr;

    user_data = (OptionDataItem *)client_data;
    
    /*
     * get values of integer parameters
     * save any changed values for log file
     */
    while (i<PARAMS_MAX && Opts_params[i]!=NULL) {
        XtVaGetValues(Opts_params[i],
               XmNuserData, &ptr,
               NULL);
        choice = *(int *)ptr;
        switch(choice) {
            case O_NX3: str = XmTextFieldGetString(Opts_params[i]);
                        val = atoi(str);
                        XtFree(str);
                        if (val!=Settings.plot_opts.nx3) {
                          sprintf(buf,"%4d",Settings.plot_opts.nx3);
                          XtVaSetValues(Opts_params[i],XmNvalue,buf,NULL);
                        }
                    break;
            case O_MP : str = XmTextFieldGetString(Opts_params[i]);
                        val = atoi(str);
                        XtFree(str);
                        if (val!=Settings.plot_opts.mp) {
                          sprintf(buf,"%4d",Settings.plot_opts.mp);
                          XtVaSetValues(Opts_params[i],XmNvalue,buf,NULL);
                        }
                    break;
            case O_MPE: str = XmTextFieldGetString(Opts_params[i]);
                        val = atoi(str);
                        XtFree(str);
                        if (val!=Settings.plot_opts.mpe) {
                          sprintf(buf,"%4d",Settings.plot_opts.mpe);
                          XtVaSetValues(Opts_params[i],XmNvalue,buf,NULL);
                        }
                    break;
            case O_NP : str = XmTextFieldGetString(Opts_params[i]);
                        val = atoi(str);
                        XtFree(str);
                        if (val!=Settings.plot_opts.np) {
                          sprintf(buf,"%4d",Settings.plot_opts.np);
                          XtVaSetValues(Opts_params[i],XmNvalue,buf,NULL);
                        }
                    break;
            case O_PRFLPTS : str = XmTextFieldGetString(Opts_params[i]);
                        val = atoi(str);
                        XtFree(str);
                        if (val!=Settings.plot_opts.profile_pts) {
                          sprintf(buf,"%4d",Settings.plot_opts.profile_pts);
                          XtVaSetValues(Opts_params[i],XmNvalue,buf,NULL);
                        }
                    break;
            case O_STIP : str = XmTextFieldGetString(Opts_params[i]);
                        val = atoi(str);
                        XtFree(str);
                        if (val!=Settings.plot_opts.stipple) {
                          sprintf(buf,"%4d",Settings.plot_opts.stipple);
                          XtVaSetValues(Opts_params[i],XmNvalue,buf,NULL);
                        }
                    break;
            case O_SOLNROT : str = XmTextFieldGetString(Opts_params[i]);
                        val = atoi(str);
                        XtFree(str);
                        if (val!=Settings.plot_opts.solution_rot) {
                          sprintf(buf,"%4d",Settings.plot_opts.solution_rot);
                          XtVaSetValues(Opts_params[i],XmNvalue,buf,NULL);
                        }
                    break;
            case O_DBLE : XmToggleButtonGadgetSetState(Opts_params[i],
                              (Settings.plot_opts.dble==1), 0);
                    break;
            case O_LABEL : XmToggleButtonGadgetSetState(Opts_params[i],
                              (Settings.plot_opts.label==1), 0);
                    break;
            case O_FLIP : XmToggleButtonGadgetSetState(Opts_params[i],
                              (Settings.plot_opts.flip==1), 0);
                    break;
            case O_CLIPTOCELL  : XmToggleButtonGadgetSetState(Opts_params[i],
                              (Settings.clip_to_cell==1), 0);
                    break;
            case O_MARKCELL : XmToggleButtonGadgetSetState(Opts_params[i],
                              (Settings.mark_cell==1), 0);
                    break;
            case O_VERBOSE  : XmToggleButtonGadgetSetState(Opts_params[i],
                              (Settings.verbose==1), 0);
                    break;
#if XY
#endif
            default:
                    break;
        }
        if (!strcmp(XtName(Opts_params[i]),"Lines"))
                XmToggleButtonGadgetSetState(Opts_params[i],
                    (Settings.plot_opts.contour_plot&LINES)!=0, 0);
        if (!strcmp(XtName(Opts_params[i]),"Shading"))
                XmToggleButtonGadgetSetState(Opts_params[i],
                    (Settings.plot_opts.contour_plot&SHADE)!=0, 0);
        i++;
    }
    XtPopdown(user_data->p_u_shell);
}

/*
 * dialog toggle button callback function
 */
void toggled_int(w, client_data, call_data)
Widget w;
XtPointer client_data, call_data;
{
    int *val = (int *)client_data;

    XmToggleButtonCallbackStruct *data = 
        (XmToggleButtonCallbackStruct *) call_data;
    if (data->set) *val = 1;
    else *val = 0;
}

void toggled_uchar(w, client_data, call_data)
Widget w;
XtPointer client_data, call_data;
{
    unsigned char *val = (unsigned char *)client_data;

    XmToggleButtonCallbackStruct *data = 
        (XmToggleButtonCallbackStruct *) call_data;
    if (data->set) *val = 1;
    else *val = 0;
}

void toggled_lines(w, client_data, call_data)
Widget w;
XtPointer client_data, call_data;
{
    int *val = (int *)client_data, mask;

    XmToggleButtonCallbackStruct *data = 
        (XmToggleButtonCallbackStruct *) call_data;
    mask = ~LINES;
    if (data->set) *val = *val | LINES;
    else *val = *val & mask;
}

void toggled_shade(w, client_data, call_data)
Widget w;
XtPointer client_data, call_data;
{
    int *val = (int *)client_data, mask;

    XmToggleButtonCallbackStruct *data = 
        (XmToggleButtonCallbackStruct *) call_data;
    mask = ~SHADE;
    if (data->set) *val = *val | SHADE;
    else *val = *val & mask;
}

void bar_cb(w, client_data, call_data)
Widget w;
XtPointer client_data, call_data;
{
    int index;
    XmString label;
    char * str;

    XtVaGetValues( w,XmNlabelString,&label,NULL );
    XmStringGetLtoR(label,"",&str);
    XmStringFree(label);
    if ((index=name_match(str,colourbar_terms))!=-1)
        Tmp.plot_opts.colour_bar = index;
    XtFree(str);
}

void set_tic(w, client_data, call_data)
Widget w;
XtPointer client_data, call_data;
{
    int index;
    XmString label;
    char * str;

    XtVaGetValues( w,XmNlabelString,&label,NULL );
    XmStringGetLtoR(label,"",&str);
    XmStringFree(label);
    if ((index=name_match(str,ticmark_terms))!=-1)
        Tmp.plot_opts.ticmark = index;
    XtFree(str);
}

static void CreateActionArea(parent, actions, num_actions, wdth)
Widget parent;
PlotActionAreaItem *actions;
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
