

依据两本书实现的基本框架，已经包含了球体和简单立方体（平行于坐标轴）两种模型。在真实渲染场景中，这还远远不够，更普遍的基本模型是三角形。三角形具有计算方便，组合灵活的优秀特性。比如，使用三角形构建一个复杂的三维图形，可以将碰撞、反射等行为的计算从三维降到二维，从而大大降低运算复杂度。为了验证接口的可扩展性并且让场景更加的丰富，我们决定增加对三角模型的支持。

我们知道，所有我们定义的碰撞模型都要去继承 $$hitable$$ 类，并实现其两个虚函数：$$hit()$$ 和$$bounding\_box()$$。$$hit()$$ 函数是用来计算光线是否与模型碰撞的，如果碰撞，还需要把碰撞点、法向量等信息记录下来。$$bounding\_box()$$ 是计算模型最小的、可以包含整个模型的立方体的函数，系统计算是否碰撞时会先看是否与这个立方体碰撞，可以显著降低运算复杂度。
