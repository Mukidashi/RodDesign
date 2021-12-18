#include "mainwindow.h"
#include <QtWidgets/QApplication>
#include<iostream>
#include<stdlib.h>
// #include <windows.h>
#include <cstdio>

int main(int argc, char *argv[])
{
	// AllocConsole();
	// FILE* fp;
	// freopen_s(&fp, "CON", "w", stdout);
	// freopen_s(&fp, "CON", "r", stdin);

	QApplication a(argc, argv);

	MainWindow mainWin;
	mainWin.show();

	return a.exec();

	// FreeConsole();
}

