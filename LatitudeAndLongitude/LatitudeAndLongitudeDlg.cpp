
// LatitudeAndLongitudeDlg.cpp: 实现文件
//

#include "stdafx.h"
#include "LatitudeAndLongitude.h"
#include "LatitudeAndLongitudeDlg.h"
#include "afxdialogex.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// 用于应用程序“关于”菜单项的 CAboutDlg 对话框

class CAboutDlg : public CDialogEx
{
public:
	CAboutDlg();

// 对话框数据
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_ABOUTBOX };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支持

// 实现
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialogEx(IDD_ABOUTBOX)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialogEx)
END_MESSAGE_MAP()


// CLatitudeAndLongitudeDlg 对话框



CLatitudeAndLongitudeDlg::CLatitudeAndLongitudeDlg(CWnd* pParent /*=NULL*/)
	: CDialog(IDD_LATITUDEANDLONGITUDE_DIALOG, pParent)
	, mImg(NULL)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CLatitudeAndLongitudeDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_BAR_X, m_bar_x);
	DDX_Control(pDX, IDC_BAR_Y, m_bar_y);
	DDX_Control(pDX, IDC_BAR_R, m_bar_r);
	DDX_Control(pDX, IDC_BAR_F, m_bar_f);
	DDX_Control(pDX, IDC_BAR_FX, m_bar_fx);
	DDX_Control(pDX, IDC_BAR_FY, m_bar_fy);
	DDX_Control(pDX, IDC_BAR_K1, m_bar_k1);
	DDX_Control(pDX, IDC_BAR_K2, m_bar_k2);
	DDX_Control(pDX, IDC_BAR_P1, m_bar_p1);
	DDX_Control(pDX, IDC_BAR_P2, m_bar_p2);
}

BEGIN_MESSAGE_MAP(CLatitudeAndLongitudeDlg, CDialog)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDLOAD, &CLatitudeAndLongitudeDlg::OnBnClickedLoad)
	ON_BN_CLICKED(IDCHECKPARAM, &CLatitudeAndLongitudeDlg::OnBnClickedCheckparam)
	ON_WM_HSCROLL()
	ON_BN_CLICKED(IDC_CHECKCORRECT, &CLatitudeAndLongitudeDlg::OnBnClickedCheckcorrect)
END_MESSAGE_MAP()


// CLatitudeAndLongitudeDlg 消息处理程序

BOOL CLatitudeAndLongitudeDlg::OnInitDialog()
{
	CDialog::OnInitDialog();

	// 将“关于...”菜单项添加到系统菜单中。

	// IDM_ABOUTBOX 必须在系统命令范围内。
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		BOOL bNameValid;
		CString strAboutMenu;
		bNameValid = strAboutMenu.LoadString(IDS_ABOUTBOX);
		ASSERT(bNameValid);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// 设置此对话框的图标。  当应用程序主窗口不是对话框时，框架将自动
	//  执行此操作
	SetIcon(m_hIcon, TRUE);			// 设置大图标
	SetIcon(m_hIcon, FALSE);		// 设置小图标

	// TODO: 在此添加额外的初始化代码
	isLoad = false;
	isCorrect = false;
	showParam = false;

	return TRUE;  // 除非将焦点设置到控件，否则返回 TRUE
}

void CLatitudeAndLongitudeDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialog::OnSysCommand(nID, lParam);
	}
}

// 如果向对话框添加最小化按钮，则需要下面的代码
//  来绘制该图标。  对于使用文档/视图模型的 MFC 应用程序，
//  这将由框架自动完成。

void CLatitudeAndLongitudeDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // 用于绘制的设备上下文

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// 使图标在工作区矩形中居中
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		CRect r = new CRect(0,0,50,50);
		// 绘制图标
		dc.DrawIcon(x, y, m_hIcon);
		//CBrush myBrush;
		//myBrush.CreateSolidBrush(RGB(255, 0, 0));
		//dc.FillRect(r, &myBrush);
	}
	else
	{
		ShowView();
		CDialog::OnPaint();
		//CDialog::UpdateWindow();// 更新windows窗口，如果无这步调用，图片显示还会出现问题
		//ShowImage(TheImage, IDC_ShowImg);// 重绘图片函数
	}
}

//当用户拖动最小化窗口时系统调用此函数取得光标
//显示。
HCURSOR CLatitudeAndLongitudeDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}



void CLatitudeAndLongitudeDlg::OnBnClickedLoad()
{
	// TODO: 在此添加控件通知处理程序代码
	CFileDialog dlg(TRUE);
	char *FilePath= 0;
	CString path;
	dlg.m_ofn.lpstrFilter = _T("*.bmp");
	if (dlg.DoModal() == IDOK)
	{
		path = dlg.GetPathName();
		FilePath = path.GetBuffer(path.GetLength());  //(LPSTR)(LPCTSTR)dlg.GetPathName();
		path.ReleaseBuffer();
		if (!(mImg = clLoadImage(FilePath)))
		{
			msg.Format(_T("加载失败"));
			AfxMessageBox(msg, MB_ICONSTOP, 0);
			return;
		}
		EasyObtainParam(&param, mImg);
		//GrayImage *gray = Image2Gray(mImg);
		//param = ObtainParam(gray);
		isLoad = true;
		f = 322.0;
		BarInit();
		ShowView();
		//msg.Format(_T("x:%d, y:%d, r:%d"), param->centerX, param->centerY, param->radius);
		//AfxMessageBox(msg);
		//UpdateData(FALSE);
	}
}

void CLatitudeAndLongitudeDlg::DisPlay()
{
	CRect rect;
	CWnd *pWin = GetDlgItem(IDC_PIC1);//获取该控件的指针，就可以对该控件直接操作了
	pWin->GetClientRect(rect);//把控件的长宽、坐标等信息保存在rect里
	int width = rect.Width();//可以获取宽和高
	int height = rect.Height();

	msg.Format(_T("高度为：%d,宽度为：%d"), height, width);
	AfxMessageBox(msg, MB_YESNO | MB_ICONSTOP, 0);

	CDC *pDc = pWin->GetDC();//获取该控件的画布
							 //有了画布，下面可以自由的画图了，想怎么画就怎么话，
	pDc->Rectangle(rect);
	CBrush myBrush;
	myBrush.CreateSolidBrush(RGB(192, 0, 233));
	pDc->FillRect(rect, &myBrush);
}

/*
1 左框 （显示原图像）  2 右框 （显示矫正图像）
true 只清屏
*/
void CLatitudeAndLongitudeDlg::ShowImage(int options, bool clear)
{
	int dlgItem;
	ClImage *img = NULL;
	CDC *pDc;
	if (options == 1)
	{
		dlgItem = IDC_PIC1;
		img = mImg;
	}
	else if (options == 2)
	{
		dlgItem = IDC_PIC2;
		img = cImg;
	}

	CWnd *pWin = GetDlgItem(dlgItem);
	pDc = pWin->GetDC();
	CRect rect;
	pWin->GetClientRect(rect);
	pDc->FillSolidRect(new CRect(1, 1, rect.Width() - 1, rect.Height() - 1), 0xF0F0F0);	//清屏
	
	if (clear)
	{
		ReleaseDC(pDc);
		return;
	}

	CDC dcCompatible;
	dcCompatible.CreateCompatibleDC(pDc);
	CBitmap bitmap;
	BYTE *pData;
	pData = new BYTE[img->width * img->height * 4];
	for (int i = 0; i<img->width*img->height * 4; i += 4)
	{
		int n = (i + 4) / 4 - 1;
		pData[i] = img->imageData[n * 3];
		pData[i + 1] = img->imageData[n * 3 + 1];
		pData[i + 2] = img->imageData[n * 3 + 2];
		pData[i + 3] = 0;
	}
	bitmap.CreateBitmap(img->width, img->height, 1, 32, pData);
	dcCompatible.SelectObject(&bitmap);
	
	int x0 = rect.Width() / 2 - img->width / 2;
	int y0 = rect.Height() / 2 - img->height / 2;
	pDc->BitBlt(x0, y0, rect.Width() - 1, rect.Height() - 1, &dcCompatible,0, 0, SRCCOPY);
	ReleaseDC(pDc);
}

// 画出param的辅助线
void CLatitudeAndLongitudeDlg::DrawParam()
{
	CRect rect;
	CWnd *pWin = GetDlgItem(IDC_PIC1);//获取该控件的指针，就可以对该控件直接操作了
	pWin->GetClientRect(rect);//把控件的长宽、坐标等信息保存在rect里
	CDC *pDc = pWin->GetDC();
	CPen gPen(PS_SOLID, 1, 0x00FF00);
	CPen rPen(PS_SOLID, 1, 0x0000FF);
	int x0 = rect.Width() / 2 - mImg->width / 2;
	int y0 = rect.Height() / 2 - mImg->height / 2;
	//圆心辅助线
	pDc->SelectObject(&rPen);
	pDc->MoveTo(x0, y0 + param.cy);
	pDc->LineTo(x0 + mImg->width, y0 + param.cy);
	pDc->MoveTo(x0 + param.cx, y0);
	pDc->LineTo(x0 + param.cx, y0 + mImg->height);
	//半径辅助线
	pDc->SelectObject(&gPen);
	pDc->MoveTo(x0 + param.cx - param.radius, y0 + param.cy - param.radius);
	pDc->LineTo(x0 + param.cx - param.radius, y0 + param.cy + param.radius);
	pDc->LineTo(x0 + param.cx + param.radius, y0 + param.cy + param.radius);
	pDc->LineTo(x0 + param.cx + param.radius, y0 + param.cy - param.radius);
	pDc->LineTo(x0 + param.cx - param.radius, y0 + param.cy - param.radius);
	ReleaseDC(pDc);
}

void CLatitudeAndLongitudeDlg::ShowView()
{
	if (!isLoad)
		return;
	ShowImage(1, false);
	if (isCorrect)
		ShowImage(2, false);
	if (showParam)
		DrawParam();
}

void CLatitudeAndLongitudeDlg::OnBnClickedCheckparam()
{
	// TODO: 在此添加控件通知处理程序代码
	if (BST_CHECKED == IsDlgButtonChecked(IDCHECKPARAM))
	{
		// 勾选
		showParam = true;
	}
	else
	{
		showParam = false;
	}
	ShowView();
}

void CLatitudeAndLongitudeDlg::BarInit()
{
	param.fx = 322.53;// / 4.26;
	param.fy = 323.13;// / 4.26;
	param.k1 = -0.037;
	param.k2 = 0.019;
	param.p1 = -0.017;
	param.p2 = 0.0049;

	m_bar_x.SetScrollRange(0, mImg->width - 1);
	m_bar_x.SetScrollPos(param.cx);
	SetDlgItemInt(IDC_EDITX, param.cx);

	m_bar_y.SetScrollRange(0, mImg->height - 1);
	m_bar_y.SetScrollPos(param.cy);
	SetDlgItemInt(IDC_EDITY, param.cy);

	m_bar_r.SetScrollRange(0, mImg->height > mImg->width ? mImg->height : mImg->width);
	m_bar_r.SetScrollPos(param.radius);
	SetDlgItemInt(IDC_EDITR, param.radius);

	m_bar_f.SetScrollRange(0, 5000);
	m_bar_f.SetScrollPos(f * 10);
	msg.Format(_T("%.1f"), f);
	SetDlgItemText(IDC_EDITF, msg);

	m_bar_fx.SetScrollRange(0, 5000);
	m_bar_fx.SetScrollPos(param.fx * 10);
	msg.Format(_T("%.1f"), param.fx);
	SetDlgItemText(IDC_EDITFX, msg);

	m_bar_fy.SetScrollRange(0, 5000);
	m_bar_fy.SetScrollPos(param.fy * 10);
	msg.Format(_T("%.1f"), param.fy);
	SetDlgItemText(IDC_EDITFY, msg);

	m_bar_k1.SetScrollRange(-5000, 5000);
	m_bar_k1.SetScrollPos(param.k1 * 10000);
	msg.Format(_T("%.4f"), param.k1);
	SetDlgItemText(IDC_EDITK1, msg);

	m_bar_k2.SetScrollRange(-5000, 5000);
	m_bar_k2.SetScrollPos(param.k2 * 10000);
	msg.Format(_T("%.4f"), param.k2);
	SetDlgItemText(IDC_EDITK2, msg);

	m_bar_p1.SetScrollRange(-2500, 2500);
	m_bar_p1.SetScrollPos(param.p1 * 10000);
	msg.Format(_T("%.4f"), param.p1);
	SetDlgItemText(IDC_EDITP1, msg);

	m_bar_p2.SetScrollRange(-2500, 2500);
	m_bar_p2.SetScrollPos(param.p2 * 10000);
	msg.Format(_T("%.4f"), param.p2);
	SetDlgItemText(IDC_EDITP2, msg);


}


void CLatitudeAndLongitudeDlg::OnHScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar)
{
	// TODO: 在此添加消息处理程序代码和/或调用默认值
	if (!isLoad)
		return;
	int level = 0;	//0: 1	1:0.1	2:0.0001
	//判断是哪个滚动条的消息
	CScrollBar *pbar = NULL;
	int editID;
	double *pdata = NULL;
	CWnd *pBarX = GetDlgItem(IDC_BAR_X);
	CWnd *pBarY = GetDlgItem(IDC_BAR_Y);
	CWnd *pBarR = GetDlgItem(IDC_BAR_R);
	CWnd *pBarF = GetDlgItem(IDC_BAR_F);
	CWnd *pBarFX = GetDlgItem(IDC_BAR_FX);
	CWnd *pBarFY = GetDlgItem(IDC_BAR_FY);
	CWnd *pBarK1 = GetDlgItem(IDC_BAR_K1);
	CWnd *pBarK2 = GetDlgItem(IDC_BAR_K2);
	CWnd *pBarP1 = GetDlgItem(IDC_BAR_P1);
	CWnd *pBarP2 = GetDlgItem(IDC_BAR_P2);
	if (pBarX == pScrollBar)
	{
		pbar = &m_bar_x;
		editID = IDC_EDITX;
		pdata = &param.cx;
	}
	if (pBarY == pScrollBar)
	{
		pbar = &m_bar_y;
		editID = IDC_EDITY;
		pdata = &param.cy;
	}
	if (pBarR == pScrollBar)
	{
		pbar = &m_bar_r;
		editID = IDC_EDITR;
		pdata = &param.radius;
	}
	if (pBarF == pScrollBar)
	{
		level = 1;
		pbar = &m_bar_f;
		editID = IDC_EDITF;
		pdata = &f;
	}
	if (pBarFX == pScrollBar)
	{
		level = 1;
		pbar = &m_bar_fx;
		editID = IDC_EDITFX;
		pdata = &param.fx;
	}
	if (pBarFY == pScrollBar)
	{
		level = 1;
		pbar = &m_bar_fy;
		editID = IDC_EDITFY;
		pdata = &param.fy;
	}
	if (pBarK1 == pScrollBar)
	{
		level = 2;
		pbar = &m_bar_k1;
		editID = IDC_EDITK1;
		pdata = &param.k1;
	}
	if (pBarK2 == pScrollBar)
	{
		level = 2;
		pbar = &m_bar_k2;
		editID = IDC_EDITK2;
		pdata = &param.k2;
	}
	if (pBarP1 == pScrollBar)
	{
		level = 2;
		pbar = &m_bar_p1;
		editID = IDC_EDITP1;
		pdata = &param.p1;
	}
	if (pBarP2 == pScrollBar)
	{
		level = 2;
		pbar = &m_bar_p2;
		editID = IDC_EDITP2;
		pdata = &param.p2;
	}


	// 获取水平滚动条当前位置
	int pos = pbar->GetScrollPos();

	switch (nSBCode)
	{
		// 如果向左滚动一列，则pos减1
	case SB_LINELEFT:
		pos -= 1;
		break;
		// 如果向右滚动一列，则pos加1
	case SB_LINERIGHT:
		pos += 1;
		break;
		// 如果向左滚动一页，则pos减10
	case SB_PAGELEFT:
		pos -= 10;
		break;
		// 如果向右滚动一页，则pos加10
	case SB_PAGERIGHT:
		pos += 10;
		break;
		// 如果滚动到最左端，则pos为1
	case SB_LEFT:
		pos = 1;
		break;
		// 如果滚动到最右端，则pos为100
	case SB_RIGHT:
		pos = 100;
		break;
		// 如果拖动滚动块到指定位置，则pos赋值为nPos的值
	case SB_THUMBPOSITION:
		pos = nPos;
		break;
		// 下面的m_horiScrollbar.SetScrollPos(pos);执行时
		// 会第二次进入此函数，最终确定滚动块位置，并且会
		// 直接到default分支，所以在此处设置编辑框中显示数值
	default:
		if (level == 1)
		{
			//SetDlgItemInt(editID, pos / 100);
			msg.Format(_T("%.1f"), *pdata);
			SetDlgItemText(editID, msg);
		}
		else if (level == 2)
		{
			msg.Format(_T("%.4f"), *pdata);
			SetDlgItemText(editID, msg);
		}
		else
			SetDlgItemInt(editID, pos);
		if (isCorrect)
			Correct();
		ShowView();
		return;
	}

	// 保存数据
	if (level == 1)
		*pdata = (double)pos / 10.0;
	else if (level == 2)
		*pdata = (double)pos / 10000.0;
	else
		*pdata = pos;
	// 设置滚动块位置
	pbar->SetScrollPos(pos);

	CDialog::OnHScroll(nSBCode, nPos, pScrollBar);
}


void CLatitudeAndLongitudeDlg::OnBnClickedCheckcorrect()
{
	// TODO: 在此添加控件通知处理程序代码
	if (!isLoad)
	{
		CButton* pBtn = (CButton*)GetDlgItem(IDC_CHECKCORRECT);
		pBtn->SetCheck(0);
		return;
	}

	if (BST_CHECKED == IsDlgButtonChecked(IDC_CHECKCORRECT))
	{
		// 勾选
		Correct();
		isCorrect = true;
	}
	else
	{
		isCorrect = false;
		ShowImage(2, true);
	}
	ShowView();
}


void CLatitudeAndLongitudeDlg::Correct()
{
	//cImg = Correct1(mImg, param, f);
	
	cImg = Correct3(mImg, &param);
}
