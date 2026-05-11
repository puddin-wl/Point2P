
loadlibrary('SecondDll.dll','SecondDll.h')%加载库文件
calllib('SecondDll','saShowImageFromFilePath','123.jpg',0,0,0,1920,1080,1)%在矩形0，0，1920，1080的位置，以不拉神模式，显示123.jpg
pause(2)%暂停2秒
calllib('SecondDll','saShowWindow',0)%隐藏窗口
pause(2)%暂停2秒
calllib('SecondDll','saShowWindow',1)%显示窗口
pause(2)%暂停2秒
calllib('SecondDll','saCloseWindow')%关闭窗口
pause(2)%暂停2秒
calllib('SecondDll','saShowImageFromFolderPath','jpg_10',0,0,0,1920,1080,500,1)%在矩形区域0，0，1920，1080的位置，以不拉伸模式，500ms为间隔，顺序播放文件夹jpg_10下的图片
pause(2)%暂停2秒
calllib('SecondDll','saPauseShow')%暂停播放
pause(2)%暂停2秒
calllib('SecondDll','saResumeShow')%恢复播放
pause(4)%暂停4秒
calllib('SecondDll','saShowImageFromFolder',0,0,0,1920,1080,500,1)%弹出一个文件夹选择窗口，可以对文件夹下的图片进行播放。在矩形区域0，0，1920，1080的位置，以不拉伸模式，500ms为间隔，顺序播放所选择的文件夹下的图片
pause(6)%暂停6秒
calllib('SecondDll','saShowImageFromSelector',0,0,0,1920,1080,500,1)%弹出一个文件过滤选择窗口，通过ctrl或shift+鼠标左键多选图片，可以对选定的图片进行播放。在矩形区域0，0，1920，1080的位置，以不拉伸模式，500ms为间隔，顺序播放所选择的图片文件。
pause(6)%暂停6秒
calllib('SecondDll','saCloseWindow')%关闭普通窗口
calllib('SecondDll','Timeout_CreateWindow')%创建带延时属性的窗口
calllib('SecondDll','Timeout_ShowWindow',1)%显示带延时属性的窗口
calllib('SecondDll','Timeout_ShowImageFromFilePath','123.jpg',0,0,0,1920,1080,1,2000)%显示图片，延时2秒
calllib('SecondDll','Timeout_ShowWindow',0)%隐藏带延时属性的窗口
calllib('SecondDll','Timeout_CloseWindow')%关闭带延时属性的窗口

%下面使用带延时属性的窗口播放图片
calllib('SecondDll','Timeout_CreateWindow')%创建带延时属性的窗口
calllib('SecondDll','Timeout_ShowWindow',1)%显示带延时属性的窗口

for i=1:10     
    str = ['jpg_10/',num2str(i),'.jpg'];
    calllib('SecondDll','Timeout_ShowImageFromFilePath',path(i),0,0,0,1920,1080,1,2000)
end
calllib('SecondDll','Timeout_CloseWindow')%关闭带延时属性的窗口

