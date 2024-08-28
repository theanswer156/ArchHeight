#### 1. 常见代码风格简述:

#### 驼峰法命名法:

第一个单词首字母小写，其它每个单词的首字母大写，其他字母小写。如:
`myNumber`   `doSomething`  

#### - 下划线命名法:

字母均小写，并将不同的单词之间用下划线 _ 分隔。如:   

`my_number`  `do_something`  

#### -  帕斯卡命名法:

将每个单词的首字母都大写。它与驼峰命名法相似，不同之处在于，帕斯卡命名法中每个单词都必须大写。如：  
`MyNumber`  `DoSomeThing`

#### -  匈牙利命名法:

以数据类型作为前缀，后跟描述变量用途的单词，且后面的单词使用帕斯卡命名法。如：

`int iCount;`   `double dPrice;`  `char cFirstChar;`

#### -  全部大写字母命名法:

将整个名称都用大写字母表示，不同单词间用下划线分隔。如：
`MY_NUMBER` `DO_SOMETHING`  

## 2. 睿达MantiSoft项目组代码规范:

#### 2.1 文件名:

文件名使用帕斯卡命名法，且在代码头文件中include的时候，也要保持和文件名一致。  
因为在Linux系统下，区分文件名大小写。  
文件中有单独的一个类， 则尽可能使文件名和类名保持一致，以便于查阅。

#### 2.2 头文件

一个头文件尽可能只包含一个类。

头文件的唯一包含使用:

```cpp
#ifndef UIPARAM_H
#define UIPARAM_H
.....代码部分.....
#endif // UIPARAM_H
```

不建议使用 `#program once` ([https://en.wikipedia.org/wiki/Pragma_once](!https://en.wikipedia.org/wiki/Pragma_once) In the C and C++ programming languages， #pragma once is a non-standard but widely supported preprocessor directive designed to cause the current header file to be included only once in a single compilation.)

使用前置声明：  
头文件中，类的成员变量声明若为指针形式，则使用前置声明，在cpp中再对此头文件进行包含。  
可提高编译速度以及减少不必要的包含文件及使用前置声明来解决循环依赖问题。

头文件中类的组织结构排序大至如下：

```cpp
class ObjName
{
public:
    ObjName();
    ~ObjName();

    ObjName &operator=(const ObjName &obj)；

public:
    void doSomething();

private:
    void setName(const std::string& str);

public:
    bool m_bFlag;

private:
    std::string m_strName;

};
```

类的构造/析构，以及符号重载函数放一起，其次为公/私有函数，公/私有变量，且注意加空行间隔。

#### 2.3 代码

#### 代码命名方式:

- 类名使用帕斯卡命名法，函数名采用驼峰命名法，变量名采用 m_匈牙利命名法，常量采用全部大写字母法。

其中，匈牙利命名法，前缀部分参考如下:

```cpp
g_ 全局变量 global
c_ 常量 const
m_ c++类成员变量 member
s_ 静态变量 static
```

基本类型部分：  
**简写 = 类型**    

```cpp
n = int 
d = double/float 
l = long 
b = bool 
c = char 

str = string   

arr = 数组  
v/vec = vector 
lst = list  
```

Qt 控件部分：

```cpp
dlg = QDialog
chk = QCheckBox 
btn = QPushButton/Button  
cmb = ComboBox 
edt = LineEdit  
spb = QSpinBox 
lbl = QLabel  
act = QAction
wnd = Window 
rb = QRadioButton  
fLayout = QFormLayout
hLayout = QHBoxLayout;
vLayout = QVBoxLayout;
gLayout = QGridLayout;
```

其它命名的简写：

```cpp
p = pointer
pt = point
func = function
img = image
alg = algorithm
num = number
it/iter = iterator
cc = cancle
del = delete
opt = option
seg = segment
2 = to
L2R = LeftToRight
B2T = BottomToTop
calc = calculate
```

若其它命名简写可能不太寻常,可在变量定义后面加上单词的全拼,使代码阅读者方便理解。  

- 函数命名规范：   
  函数命名应做到“一字见心”的原则,如下示例:
  do：表示执行一个动作，例如 `doSomething()`。  
  is：表示查询某个值，例如 `isEnable()`。  
  get：表示获取某个值，例如 `getName()`。  
  set：表示设置某个值，例如 `setName()`。  

信号与槽的函数，前缀分别加上‘sig’， ‘slot’,如：

```cpp
void sigValueChange(int nValue); // 信号函数 
void slotApply();    // 槽函数
```

- 变量初始化
  在类的指针成员变量中，均需将指针初始化为nullptr值，
  初始化的方式可采用 = 及 {} 两种方式， 但同一个类中，必须使用同一种风格，如：
  
  ```cpp
  QSpinBox *m_spbXAxis = nullptr;  // 赋值初始化方式
  ```

QSpinBox *m_spbXAxis{ nullptr }; // C++11 统一初始化方式

#### 空格及换行：

- 代码不能数十行挤在一起，约十五至二十行左右，加一个换行符，使代码有段落感，且代码行不应过长。

- 代码行间，需加空格或使其更易阅读，如：
  
  ```cpp
  if(a>b) -> if (a > b)
  a+b  ->  a + b
  setValue(1,2,3)  -> setValue(1, 2, 3)
  ```

- if语句后若紧接一条语句，此语句需放至下一行，而不和if在同一行，在同一行使调试不方便。
  if与else间，若代码有些为一行，有些为若干行，则均需加{}，使代码前后观感保持一致，若均为一行，可加可不加。
  
  ```cpp
  if (bStatus)
  {
    行1
    行2
  }
  else
  {
    行1
  }
  ```
  
  for 和 if 在一起时,for的语句部分必须加上{},如:
  
  ```cpp
  // 难以阅读的写法
  for(const auto& item : vItems)
    if(item)
        item->setName("Test");
  ```  
  
```cpp
// 更为美观的写法
for(const auto& item : vItems)
{
    if(item)
        item->setName("Test");
}

```
- 函数与函数间，需添加空行。

- 大括号均放至一条竖线上

- 代码尽可能使用double，而不使用float类型，以避免精度损失问题。

- 不允许采用拼音方式对变量进行命名。
- 多个同类型的变量，不允许使用数字进行标识。
- 不允许使用 MyXXX， 推荐使用 RdXXX。

- 编写界面参数，在类中，变量按照界面的控件由上至下进行排列，以方便参照界面进行函数/变量的查找。

- 代码编写需整洁美观，有段落感。有若干行时，如许多变量/函数调用间，加空格，以使代码错落有致。

- 推荐VS编译器设置如下：
在VS菜单中，依次到“工具->选项->文本编辑器->C/C++->代码样式->格式设置->新行->”，进行设置：
1. 命名空间的左大括号的位置，类型的左大括号的位置，函数的左大括号的位置，控制块的左大括号的位置，均设置为“移动到新行”。  
2. lambda的左括号的位置 设置为 保持在同一行上，但在前面添加一个空格。  
3. 其它：  
3.1 将"catch"和相似的关键字放在新行上。  
3.2 将"else"放在新行上。  
使用：按住Ctrl 然后依次按下 K，F键， 可对选择的代码进行格式化，或在菜单中“编辑->高级->设置选定内容的格式"，对代码进行格式化。

亦或在“扩展”菜单中，安装“CodeMaid”控件，插件使用VS的设置格式，对代码批量进行格式化，并删除多余空行。  

#### 其它约束：   
- 信号槽连接，使用Qt5的连接方式:
```cpp
// Qt4的连接方式：
QObject::connect(sender, SIGNAL(signal()), receiver, SLOT(slot()));

// Qt5的连接方式：
QObject::connect(sender, &Sender::signal, receiver, &Receiver::slot);

// Qt5的连接方式（有重载的信号或槽）：
QObject::connect(btn, &QPushButton::clicked, 
    receiver, static_cast<void(RdClass::*)(int)>(&RdClass::slotBtnClick));
```

若界面类,控件的成员变量及操作函数,应按界面中的视觉排序方式从左到右从上到下,依次编写.使用者查询界面中的元素更为方便.

编译器编译过程中的警告提示,尽可能的进行解决,避免意外的bug
声明指针必须初始化为空.


#### 2.4 注释

- 对一个类，应简略描述其作用。算法类或数据结构,应详细描述算法的大致过程，及数据的使用场景。

- 若函数比较复杂,应写上详细的注释说明, 一些简单的函数注释可使用一至两行写些简单的注释， 如:
  ```double getMaxNum();    // 获取最大数值```
  注释的//符和代码间，建议加tab或空格，使其间隔开来。

- 推荐头文件的注释方式：
  
  ```cpp
  /*
  * @Author: 
  * @Date: 
  * @LastEditors: 
  * @LastEditTime:
  * @Description: 
  * 
  */
  ```

- 推荐函数的注释方式：
  
  ```cpp
  /***
  * @brief        : none
  * @in_param : none
  * @out_param : none
  * @return : none
  */
  ```
  
  - 代码标签，代码中特殊的注释技术
    代码中,部分函数可能未完成或有bug等情况,可添加一些特殊的标注,进行备忘。
    标签使用说明如下：
    
    ```cpp
    // TODO：标记代码中需要实现的功能或任务。
    // FIXME：标记代码中需要修复的问题或缺陷。
    // NOTE：提供额外的注释或提示信息，帮助理解代码意图或设计决策。
    // BUG：标记已知的Bug或错误。
    // XXX：标记需要警惕或需要重点关注的代码块。
    // HACK：标记临时性修复或不优雅的解决方案。
    ```
    
    在VS中，可通过“视图->其他窗口->任务列表”，将上述标签显示出来。
    VSCode中，可安装“Todo Tree”等插件，显示这些特殊的标签。
 
 代码中尽量少用 /* */ 对一大串代码进行注释，因为 /*这种注释嵌套使用容易出问题，
 可使用 #if 0 #endif 代替

  #### 2.5 Git提交规范
  
  Git常见的提交类型包括：  
  
  ```cpp
  feat：新功能
  fix：修复 bug
  docs：文档修改
  style：代码格式修改，比如缩进、空格等
  refactor：代码重构
  test：测试相关修改
  chore：其他修改，比如构建流程、辅助工具等  
  perf: 优化相关，比如提升性能、体验
  revert: 回滚到上一个版本
  ```

推荐使用: 

```cpp
新增 --- 添加了计算曲线曲率的函数接口
修复 --- 图元拖拽丢失
重构 --- 重写了计算面积的函数接口
调整 --- 对XXX类重命名及代码调整
优化 --- 提升了XXX算法的效率
回退 --- 因前版本XXX的Bug，回退至XXX并重写

即 大概功能 间隔符  大体描述, 可自定义编写
```

若修改较大，且描述较多，则使用多行注释，注释大体格式如下：

```cpp
# 标题行 类型（影响范围）：50个字符以内，描述主要变更内容
# 空一行
# 主体内容：更详细的说明文本，建议72个字符以内。 需要描述的信息包括:
# * 为什么这个变更是必须的? 它可能是用来修复一个bug，增加一个feature，提升性能、可靠性、稳定性等等
# * 他如何解决这个问题? 具体描述解决问题的步骤
# * 是否存在副作用、是否会有风险? 
# 空一行
# 尾部：是需要的化可以添加一个链接到issue地址或者其它文档，或者是关闭了某个issue。
```

示例:

```cpp
新增  --- 添加XXX算法

算法描述1
算法描述2
算法描述3
```  
提交尽量要功能单一,不应在一次提交中包含大量功能的修改



#### 代码优化

使用const引用传递而非值传递通过值传递给函数时会创建临时变量。
`void process(vector<int> s);`  
改为引用传递：  
`void process(const vector<int>& s);`  
如果参数对象较大，使用值传递时拷贝的代价大，应该使用引用传递。  
如果是内置类型int, float, double, char 和 bool，使用值传递，可以支持某一些编译器优化手段。  
同理，对类成员访问接口：
```cpp
Status GetStates() const
{
    return m_status;
}
```
应改为返回 const&，
```cpp
const Status& GetStates() const
{
    return m_status;
}
```
保存函数返回的引用时，为避免拷贝，同样需要引用：`const auto& states = GetStates();`

返回引用避免拷贝还需要避免隐式转换的发生，否则还是会拷贝。

for循环中使用引用遍历以下for循环，每次从 object_list 容器拷贝临时对象 obj,

```cpp
for(auto obj: object_list) 
{
 obj->func();
}
```
改为 const auto& 避免拷贝。
```cpp
for(const auto& obj: object_list) 
{
 obj->func();
}
```