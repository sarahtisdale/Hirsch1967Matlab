function varargout = HirschGraphs(varargin)
% HIRSCHGRAPHS MATLAB code for HirschGraphs.fig
%      HIRSCHGRAPHS, by itself, creates a new HIRSCHGRAPHS or raises the existing
%      singleton*.
%
%      H = HIRSCHGRAPHS returns the handle to a new HIRSCHGRAPHS or the handle to
%      the existing singleton*.
%
%      HIRSCHGRAPHS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HIRSCHGRAPHS.M with the given input arguments.
%
%      HIRSCHGRAPHS('Property','Value',...) creates a new HIRSCHGRAPHS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HirschGraphs_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HirschGraphs_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HirschGraphs

% Last Modified by GUIDE v2.5 08-Feb-2016 20:13:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @HirschGraphs_OpeningFcn, ...
                   'gui_OutputFcn',  @HirschGraphs_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before HirschGraphs is made visible.
function HirschGraphs_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HirschGraphs (see VARARGIN)

% Choose default command line output for HirschGraphs
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Update inner K+/lambda+ controls
CheckInnerAutoChange(handles);
first_plots(hObject, handles);

% UIWAIT makes HirschGraphs wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function first_plots(hObject, handles)
% first plot is sets up initial plot so all other can be "replots" 

%%%% fig3 - initial plot and one-time attributes
axes(handles.AxesFig3);
handles.fig3plot = plot([1 2],[1 1]); % dummy plot for now, real data in fig3plot()
handles.fig3text = text(0,0,'foo');   % dummy label for now, real label fig3plot()
guidata(hObject, handles);
ylabel('\phi(R)');
xlabel('R');
fig3_new_radius(handles);

%%%% fig5 - initial plot and one-time attributes
axes(handles.AxesFig5);
handles.fig5plotI = plot([0 1],[1 1]); % dummy plot for now, real data in fig5plot()
handles.fig5plotO = plot([1 2],[1 1]); % dummy plot for now, real data in fig5plot()
guidata(hObject, handles);
ylabel('\phi(R)');
xlabel('R');

%%%% fig6 - initial plot and one-time attributes
axes(handles.AxesFig6);
handles.fig6plot = plot([1 2],[1 1]); % dummy plot for now, real data in fig6plot()
handles.fig6text = text(0,0,'foo');   % dummy label for now, real label fig6plot()
guidata(hObject, handles);
ylabel('Normalized Ion Density');
xlabel('R');

%%%% fig7 - initial plot and one-time attributes
axes(handles.AxesFig7);
handles.fig7plotP = plot([0 1],[1 1]); % dummy plot for now, real data in fig7plot()
R = 0:0.005:1;
plot(R,0.1122./R.^2,'--');
handles.fig7plotE = plot([1 2],[1 1]); % dummy plot for now, real data in fig7plot()
guidata(hObject, handles);
ylabel('Normalized Pressure');
xlabel('R');
legend_strs = cell(3,1);
legend_strs{1} = 'P / V_{0}^{2}';
legend_strs{2} = '0.1122 / R^{2}';
legend_strs{3} = 'E^{2} / 8\piV_{0}^{2}';
legend(legend_strs,'Location','southeast'); % add legend

replot_figs(handles);

function replot_figs(handles)
fig3plot(handles);
fig567plots(handles);

function fig3_new_radius(handles)
%%%% fig 3
axes(handles.AxesFig3);
ax = gca;
RMax = str2double(handles.EditFig3R.String);
ax.XTick = 0:1:RMax;
ax.XLim = [0 RMax];


function fig3plot(handles)
% radius from 1 to 4
Rmax = str2double(handles.EditFig3R.String);
R = 1:0.01:Rmax;

% intial values at radius 1 (from the text):
hinits=[0.0000001,0];

K_plus      = str2double(handles.EditOuterK.String);
lambda_plus = str2double(handles.EditOuterL.String);

% solving
[R,v] = ode45(@hirsch_equ7_laplace,R,hinits,[],[K_plus lambda_plus]);

% strip imaginary values (can't plot them)
v(find(real(v)~=v)) = NaN;

% re-plot
set(handles.fig3plot,'XData',R,'YData',v(:,1));  % only plot the first column (second column is v')

% manually find maximum
maxv = max(v(:,1));
maxidx = find(v(:,1)==maxv);
rmax = R(maxidx);
% create a silly text label
% alignment='left';
%str = sprintf('R_{vmax}=%f, v=%f',rmax,maxv);
%label_str = sprintf('\\leftarrow %s',str);
% if maxidx > (numel(R)/2)
%     alignment='right';
%    label_str = sprintf('%s \\rightarrow',str);
% end
%set(handles.fig3text,'String',label_str,'Position',[R(maxidx) maxv],'HorizontalAlignment',alignment);

% Show Rc/Ra in text field
handles.RcOverRa.String = num2str(rmax);
handles.RaOverRc.String = num2str(1/rmax);

function fig567plots(handles)

%%%%%%%%%%%%
% fig 5

% outer K+=0.7 plot
%  - solved backward from real cathode to virtual anode 
K_plus      = str2double(handles.EditOuterK.String);
lambda_plus = str2double(handles.EditOuterL.String);
% radius varies from 1 down to virtual anode
Ra_over_Rc = str2double(handles.RaOverRc.String);
R_c = 1:-0.005:Ra_over_Rc;
% initial values: potential 1, slope zero
hinits=[0.999999999;0]; 
% solve
[R,v] = ode45(@hirsch_equ7_laplace,R_c,hinits,[],[K_plus lambda_plus]);
% strip imaginary values (can't plot them)
v(find(real(v)~=v)) = NaN;
% re-plot
set(handles.fig5plotO,'XData',R,'YData',v(:,1));  % only plot the first column (second column is v')
Rsave1 = R;
Vsave1 = v;

% inner K+=0.0859 plot
%  - solved backward from virtual anode to zero
K_plus      = str2double(handles.EditInnerK.String);
lambda_plus = str2double(handles.EditInnerL.String);
% radius varies from virtual anode down to zero;
R_c = Ra_over_Rc:-0.005:0;
% initial values : potenial and slope zero
hinits=[0.0000001;0]; 
% solve
[R,v] = ode45(@hirsch_equ7_laplace,R_c,hinits,[],[K_plus lambda_plus]);
% strip imaginary values (can't plot them)
v(find(real(v)~=v)) = NaN;
% re-plot
set(handles.fig5plotI,'XData',R,'YData',v(:,1));  % only plot the first column (second column is v')
Rsave2=R;
Vsave2=v;


%%% figure 6

% set up full R&V as from two previous plots
R = [flipud(Rsave2)' flipud(Rsave1)'];
V = [flipud(Vsave2(:,1))' flipud(Vsave1(:,1))'];
E = [flipud(Vsave2(:,2))' flipud(Vsave1(:,2))'];

% caluclate ion density
%  - use equ 8 for K+, but re-solve for rho_i
V_0 = 1; % V_0 in our normalized world is 1
K_plus = 0.7;
rho_i = abs(V_0).*K_plus./(4.*pi.*(R.^2).*(V.^(1/2)));
set(handles.fig6plot,'XData',R,'YData',rho_i);

%%%%%%%%%%%%%
% fig 7

% set up graph, axis labels and limits

% caluclate ion pressure
%  - use equ 17, but note the erratum! 
V_0 = 1; % V_0 in our normalized world is 1
r_c = 1;
P = ((K_plus.*V_0.^2.)/(2.*pi.*r_c.^2)) .* (1./R.^2) .* ...
    (lambda_plus.*((1-V).^(1./2))+V.^(1./2));
set(handles.fig7plotP,'XData',R,'YData',P);
set(handles.fig7plotE,'XData',R,'YData',E.^2/(8.*pi));





% --- Outputs from this function are returned to the command line.
function varargout = HirschGraphs_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function EditOuterK_Callback(hObject, eventdata, handles)
% hObject    handle to EditOuterK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Log10EditCB(hObject,handles.SliderOuterK);
AutoInnerUpdate(handles);
replot_figs(handles);

% --- Executes during object creation, after setting all properties.
function EditOuterK_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditOuterK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function SliderOuterK_Callback(hObject, eventdata, handles)
% hObject    handle to SliderOuterK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
newK = get(hObject,'Value');
newK = 10^newK;
set(handles.EditOuterK, 'String', num2str(newK,4));
AutoInnerUpdate(handles);
replot_figs(handles);


% --- Executes during object creation, after setting all properties.
function SliderOuterK_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SliderOuterK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function EditOuterL_Callback(hObject, eventdata, handles)
% hObject    handle to EditOuterL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
EditCB(hObject,handles.SliderOuterL);
AutoInnerUpdate(handles);
replot_figs(handles);


% --- Executes during object creation, after setting all properties.
function EditOuterL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditOuterL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function SliderOuterL_Callback(hObject, eventdata, handles)
% hObject    handle to SliderOuterL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
newL = get(hObject,'Value');
set(handles.EditOuterL, 'String', num2str(newL,4));
AutoInnerUpdate(handles);
replot_figs(handles);


% --- Executes during object creation, after setting all properties.
function SliderOuterL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SliderOuterL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function EditInnerK_Callback(hObject, eventdata, handles)
% hObject    handle to EditInnerK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Log10EditCB(hObject,handles.SliderInnerK);
replot_figs(handles);



% --- Executes during object creation, after setting all properties.
function EditInnerK_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditInnerK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function SliderInnerK_Callback(hObject, eventdata, handles)
% hObject    handle to SliderInnerK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
newK = get(hObject,'Value');
newK = 10^newK;
set(handles.EditInnerK, 'String', num2str(newK,4));
replot_figs(handles);


% --- Executes during object creation, after setting all properties.
function SliderInnerK_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SliderInnerK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function EditInnerL_Callback(hObject, eventdata, handles)
% hObject    handle to EditInnerL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
EditCB(hObject,handles.SliderInnerL);
replot_figs(handles);

% --- Executes during object creation, after setting all properties.
function EditInnerL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditInnerL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function SliderInnerL_Callback(hObject, eventdata, handles)
% hObject    handle to SliderInnerL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
newL = get(hObject,'Value');
set(handles.EditInnerL, 'String', num2str(newL,4));
replot_figs(handles);

% --- Executes during object creation, after setting all properties.
function SliderInnerL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SliderInnerL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function Log10EditCB(edit,slider)
new = edit.String;
new = str2double(new);
new_edit = new;
new = log10(new);
if (isnan(new) || new < slider.Min || new > slider.Max)
    new_edit = 10^slider.Value;
    new = slider.Value;
    errordlg(sprintf('Enter a number between %.4f and %.4f.',10^slider.Min,10^slider.Max));
end
set(edit, 'String', num2str(new_edit));
set(slider, 'Value', new);


function EditCB(edit,slider)
new = edit.String;
new = str2double(new);
if (isnan(new) || new < slider.Min || new > slider.Max)
    new = slider.Value;
    errordlg(sprintf('Enter a number between %.4f and %.4f.',slider.Min,slider.Max));
end
set(edit, 'String', num2str(new));
set(slider, 'Value', new);


function EditFig3R_Callback(hObject, eventdata, handles)
% hObject    handle to EditFig3R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
EditCB(hObject,handles.SliderFig3R);
fig3_new_radius(handles);
fig3plot(handles);

% --- Executes during object creation, after setting all properties.
function EditFig3R_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditFig3R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function SliderFig3R_Callback(hObject, eventdata, handles)
% hObject    handle to SliderFig3R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
newR = get(hObject,'Value');
set(handles.EditFig3R, 'String', num2str(newR));
fig3plot(handles);
fig3_new_radius(handles);


% --- Executes during object creation, after setting all properties.
function SliderFig3R_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SliderFig3R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function AutoInnerUpdate(handles)
%  Values calculated by comparing inner/outer K+/lambda+ values
%  from fig5 in Hirsch 1967 ERRATUM.  Only one set, so guessing a
%  proportional scaling, but could be (probably is) more complex.
if (~handles.CheckInnerAuto.Value)
    return;
end
Kfactor = 0.0859/0.7;
Lfactor = 3.7/0.454;
set(handles.EditInnerK, 'String', Kfactor*str2double(handles.EditOuterK.String));
set(handles.EditInnerL, 'String', Lfactor*str2double(handles.EditOuterL.String));
set(handles.SliderInnerK, 'Value', log10(Kfactor*str2double(handles.EditOuterK.String)));
set(handles.SliderInnerL, 'Value', log10(Lfactor*str2double(handles.EditOuterL.String)));


function CheckInnerAutoChange(handles)
% if enabled, calc inner K+/lambda+ values and disable edits/sliders
auto = get(handles.CheckInnerAuto,'Value');
enable_str = 'on';
if (auto)
    enable_str = 'off';
    % calc inner k+/lambda+
    AutoInnerUpdate(handles);
end
% enable/disable inner K+/lambda+ edits/sliders
set(handles.EditInnerK, 'Enable', enable_str);
set(handles.EditInnerL, 'Enable', enable_str);
set(handles.SliderInnerK, 'Enable', enable_str);
set(handles.SliderInnerL, 'Enable', enable_str);


% --- Executes on button press in CheckInnerAuto.
function CheckInnerAuto_Callback(hObject, eventdata, handles)
CheckInnerAutoChange(handles);
replot_figs(handles);

function AxesMaxMin(axes,zoompanel)
global saved_axes_parent_panel;
panel = get(axes,'Parent');
if (panel == zoompanel)
    set(axes,'parent',saved_axes_parent_panel);
    set(zoompanel,'Visible','off');
else
    saved_axes_parent_panel = panel;
    set(axes,'parent',zoompanel);
    set(zoompanel,'Visible','on');
    set(zoompanel,'Title',panel.Title);
end


% --- Executes on mouse press over axes background.
function AxesButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to AxesFig3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AxesMaxMin(hObject,handles.zoompanel);


% --- Executes during object creation, after setting all properties.
function AxesFig3_CreateFcn(hObject, eventdata, handles)
