
#ifdef _WIN32

#include <Rcpp.h>

#include <string>
#include <sstream>
#include <stdexcept>
using namespace std;

#include "wingui.h"


class findR{
  private:
    DWORD _Rpid;
    HWND  _RWnd;
  public:
    findR(DWORD pid):_Rpid(pid),_RWnd(NULL){}
    bool isR( HWND in) {
        DWORD inPid=0;
        GetWindowThreadProcessId(in, &inPid);
        // Rprintf("ipPid=%d\t", inPid);
        if (inPid == _Rpid) {
            _RWnd = in;
            // Rprintf("this is R(%p).\n", _RWnd);
            return (inPid == _Rpid);
        }
        return false;
    }
    HWND RWnd(){
        // Rprintf("RWnd=%p\t", _RWnd);
        return _RWnd;
    }
};
BOOL CALLBACK EnumFind(HWND aWnd, LPARAM lParam) {
	findR* find = (findR *)lParam;
	if (!IsWindowVisible(aWnd)){ // Skip hidden windows.
        // Rprintf("Awnd=%p\n", aWnd);
		return true;
    }
	return !(find->isR(aWnd));
}

int WindowsGUI::get_pid() {
    return GetCurrentProcessId();
}
WindowsGUI::WindowsGUI() {
    findR find(get_pid());
    EnumWindows(EnumFind, (LPARAM)&find);
    HGUI = find.RWnd();
}
string WindowsGUI::get_win() {
    ostringstream stringStream;
    stringStream << HGUI;
    return stringStream.str();
}
string WindowsGUI::get_window_text() {
    const int nMaxCount = 255;
    char szBuffer[nMaxCount+1];
    GetWindowText( HGUI, szBuffer, nMaxCount );
    string title = szBuffer;
    return title;
}
void WindowsGUI::set_window_text(string title) {
    SetWindowText(HGUI, title.c_str());
    return;
}
bool WindowsGUI::get_on_top(){
    return (GetWindowLong(HGUI, GWL_EXSTYLE) & WS_EX_TOPMOST) == WS_EX_TOPMOST;
}
void WindowsGUI::set_on_top(bool value){
    if(value){
        SetWindowPos( HGUI, HWND_TOPMOST
                    , 0,0,0,0
                    , SWP_NOREPOSITION | SWP_NOSIZE | SWP_NOMOVE
                    );
    } else {
        SetWindowPos( HGUI, HWND_NOTOPMOST
                    , 0,0,0,0
                    , SWP_NOREPOSITION | SWP_NOSIZE | SWP_NOMOVE
                    );
    }
    
}

bool WindowsGUI::get_layered(){
    return GetWindowLong(HGUI, GWL_EXSTYLE) & WS_EX_LAYERED;
}
void WindowsGUI::set_layered(bool value){
    LONG exstyle=GetWindowLong(HGUI, GWL_EXSTYLE);
    if(value)
        SetWindowLong(HGUI, GWL_EXSTYLE, exstyle | WS_EX_LAYERED);
    else
        SetWindowLong(HGUI, GWL_EXSTYLE, exstyle &~WS_EX_LAYERED);
}


double WindowsGUI::get_opacity(){
    if(GetWindowLong(HGUI, GWL_EXSTYLE) & WS_EX_LAYERED){
		COLORREF color;
		BYTE alpha;
		DWORD flags;
        GetLayeredWindowAttributes(HGUI, &color, &alpha, &flags);
        return alpha / 255.;
    } else {
        return 1;
    }
}
void WindowsGUI::set_opacity(double percent){
    if (percent < 0) percent = 0;
    else {
        LONG exstyle = GetWindowLong(HGUI, GWL_EXSTYLE);
        if (percent >= 1) {
            SetWindowLong(HGUI, GWL_EXSTYLE, exstyle & ~WS_EX_LAYERED);
        } else {
            SetWindowLong(HGUI, GWL_EXSTYLE, exstyle | WS_EX_LAYERED);
            SetLayeredWindowAttributes(HGUI, 0, 255 * percent, LWA_ALPHA);
        }
    }
}

double WindowsGUI::get_transparency(){
    return 1-get_opacity();
}
void WindowsGUI::set_transparency(double percent){
    if(percent < 0) percent = 0;
    if(percent > 1) percent = 1;
    set_opacity(1-percent);
}

#endif
