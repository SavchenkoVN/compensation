# Compensation project
This code tries to solve the problem of error compensation in KTRV

These are my observations
- The best method is that is used by @scipy fsolve (Ludwig - smthing bla bla)
- Need to think about starting point some time it's better to start with factical S21, some time need to start from modeling data

## Math section
Уравнение биений:
$$
\begin{equation}
    \begin{cases}
    A_{rx}^{2} = A_{src}^{2} + A_{inc}^{2} + 2 A_{src}A_{inc}\cos{(\varphi_{src} - \varphi_{inc})}\\
    \tan{\varphi_{rx}} = \frac{A_{src}\sin{\varphi_{src}} + A_{inc}\sin{\varphi_{inc}}}{A_{src}\cos{\varphi_{src}} + A_{inc}\cos{\varphi_{inc}}}
\end{cases}
\end{equation}
$$

Из-за плохой синхронизации SDR, а также отсутствия синхронизации между приемником и передатчиком необходимо переходить от абсолютной фазы к относительной

$$
\begin{equation}
    \begin{pmatrix} \Delta_{rx1} \\ \Delta_{rx2} \\ \Delta_{rx3} \\ \Delta_{rx4}\end{pmatrix} = \begin{pmatrix}
\varphi_{A2} \\ \varphi_{A1} \\ \varphi_{A0} \\ \varphi_{A3}
\end{pmatrix} -\begin{pmatrix}
\varphi_{A3} \\ \varphi_{A2} \\ \varphi_{A2} \\ \varphi_{A0}
\end{pmatrix} = 
\begin{pmatrix}
    \varphi_{A2} - \varphi_{A2.p0} \\
    \varphi_{A1} - \varphi_{A1.p0} \\
    \varphi_{A0} - \varphi_{A0.p0} \\
    \varphi_{A3} - \varphi_{A3.p0}
\end{pmatrix} - 
\begin{pmatrix}
    \varphi_{A3} - \varphi_{A3.p0} \\
    \varphi_{A2} - \varphi_{A2.p0} \\
    \varphi_{A1} - \varphi_{A1.p0} \\
    \varphi_{A0} - \varphi_{A0.p0}
\end{pmatrix} = 
\end{equation}

$$
$$
=
\begin{pmatrix}
    \varphi_{A2} - \varphi_{A2.p0} - \varphi_{A3} + \varphi_{A3.p0} \\
    \varphi_{A1} - \varphi_{A1.p0} - \varphi_{A2} + \varphi_{A2.p0} \\
    \varphi_{A0} - \varphi_{A0.p0} -  \varphi_{A1} + \varphi_{A1.p0} \\
    \varphi_{A3} - \varphi_{A3.p0} - \varphi_{A0} + \varphi_{A0.p0}
\end{pmatrix} = 
\begin{pmatrix}
    \varphi_{A2} - \varphi_{A3} \\
    \varphi_{A1} - \varphi_{A2} \\
    \varphi_{A0} - \varphi_{A1} \\
    \varphi_{A3} - \varphi_{A0} 
\end{pmatrix} + 
\begin{pmatrix}
    \varphi_{A3.p0} - \varphi_{A2.p0} \\
    \varphi_{A2.p0} - \varphi_{A1.p0} \\
    \varphi_{A1.p0} - \varphi_{A0.p0} \\
    \varphi_{A0.p0} - \varphi_{A3.p0} 
\end{pmatrix}
$$
Раскрываем понятия фазы через мнимую и действительную часть. Принимаем, что 
$$
\begin{equation}
    \begin{pmatrix}
    \varphi_{A3.p0} - \varphi_{A2.p0} \\
    \varphi_{A2.p0} - \varphi_{A1.p0} \\
    \varphi_{A1.p0} - \varphi_{A0.p0} \\
    \varphi_{A0.p0} - \varphi_{A3.p0} 
\end{pmatrix} = error_{const}
\end{equation}
$$
Тогда 
$$
\begin{pmatrix}
    \varphi_{A2} - \varphi_{A3} \\
    \varphi_{A1} - \varphi_{A2} \\
    \varphi_{A0} - \varphi_{A1} \\
    \varphi_{A3} - \varphi_{A0} 
\end{pmatrix} = 
\begin{pmatrix}
    angle(I_{A2} + jQ_{A2}) - angle(I_{A3} + jQ_{A3}) \\
    angle(I_{A1} + jQ_{A2}) - angle(I_{A2} + jQ_{A2}) \\
    angle(I_{A0} + jQ_{A0}) - angle(I_{A1} + jQ_{A1}) \\
    angle(I_{A3} + jQ_{A4}) - angle(I_{A0} + jQ_{A0}) 
\end{pmatrix}
$$
Раскрываем, через
$$
\begin{equation}
    S_{Ai} = I_{Ai} + jQ_{Ai} = (I_{SRC} + I_{Incident}) + j(Q_{SRC} + Q_{Incident})
\end{equation}
$$
В результате получаем 
$$
\begin{equation}
\begin{pmatrix} \Delta_{rx1} \\ \Delta_{rx2} \\ \Delta_{rx3} \\ \Delta_{rx4}\end{pmatrix} =
\begin{pmatrix}
    angle(I_{SRC.A3} + I_{Incident.A3} + j(Q_{SRC.A3} + Q_{Incident.A3})) - angle(I_{SRC.A2} + I_{Incident.A2} + j(Q_{SRC.A2} + Q_{Incident.A2})) \\
    angle(I_{SRC.A2} + I_{Incident.A2} + j(Q_{SRC.A2} + Q_{Incident.A2})) - angle(I_{SRC.A1} + I_{Incident.A1} + j(Q_{SRC.A1} + Q_{Incident.A1})) \\
    angle(I_{SRC.A1} + I_{Incident.A1} + j(Q_{SRC.A1} + Q_{Incident.A1})) - angle(I_{SRC.A0} + I_{Incident.A0} + j(Q_{SRC.A0} + Q_{Incident.A0})) \\
    angle(I_{SRC.A0} + I_{Incident.A0} + j(Q_{SRC.A0} + Q_{Incident.A0})) - angle(I_{SRC.A3} + I_{Incident.A3} + j(Q_{SRC.A3} + Q_{Incident.A3})) \\
\end{pmatrix} + error_{const}
\end{equation}
$$
4 уравнения - 8 неизвестных. Дополним амплитудой
$$
\begin{equation}
\begin{matrix}
    A_{A0}^2 = (I_{SRC.A0} + I_{Incident.A0})^2 + (Q_{SRC.A0} + Q_{Incident.A0})^2\\
    A_{A1}^2 = (I_{SRC.A1} + I_{Incident.A1})^2 + (Q_{SRC.A1} + Q_{Incident.A1})^2\\
    A_{A2}^2 = (I_{SRC.A2} + I_{Incident.A2})^2 + (Q_{SRC.A2} + Q_{Incident.A2})^2\\
    A_{A3}^2 = (I_{SRC.A3} + I_{Incident.A3})^2 + (Q_{SRC.A3} + Q_{Incident.A3})^2
\end{matrix}
\end{equation}
$$
В результате получаем 8 уравнений и 8 неизвестных
$$
\begin{equation}
\begin{cases}
    \varphi_{A2} - \varphi_{A3} =  angle(I_{SRC.A3} + I_{Incident.A3} + j(Q_{SRC.A3} + Q_{Incident.A3})) - angle(I_{SRC.A2} + I_{Incident.A2} + j(Q_{SRC.A2} + Q_{Incident.A2}))   \\
    \varphi_{A1} - \varphi_{A2} = angle(I_{SRC.A2} + I_{Incident.A2} + j(Q_{SRC.A2} + Q_{Incident.A2})) - angle(I_{SRC.A1} + I_{Incident.A1} + j(Q_{SRC.A1} + Q_{Incident.A1}))    \\
    \varphi_{A0} - \varphi_{A1} =  angle(I_{SRC.A1} + I_{Incident.A1} + j(Q_{SRC.A1} + Q_{Incident.A1})) - angle(I_{SRC.A0} + I_{Incident.A0} + j(Q_{SRC.A0} + Q_{Incident.A0}))    \\
    \varphi_{A3} - \varphi_{A0} = angle(I_{SRC.A0} + I_{Incident.A0} + j(Q_{SRC.A0} + Q_{Incident.A0})) - angle(I_{SRC.A3} + I_{Incident.A3} + j(Q_{SRC.A3} + Q_{Incident.A3}))  \\
    A_{A0}^2 = (I_{SRC.A0} + I_{Incident.A0})^2 + (Q_{SRC.A0} + Q_{Incident.A0})^2\\
    A_{A1}^2 = (I_{SRC.A1} + I_{Incident.A1})^2 + (Q_{SRC.A1} + Q_{Incident.A1})^2\\
    A_{A2}^2 = (I_{SRC.A2} + I_{Incident.A2})^2 + (Q_{SRC.A2} + Q_{Incident.A2})^2\\
    A_{A3}^2 = (I_{SRC.A3} + I_{Incident.A3})^2 + (Q_{SRC.A3} + Q_{Incident.A3})^2
\end{cases} + error
\end{equation}
$$
На основе уравнения (7) мы можем найти помеху для каждого канала.
Получаем значения амплитуды и фазы помехи в пространстве сканера. Перейдем к амлпитудам и фазам, согласно формуле:
$$
\begin{equation}
    \begin{cases}
        A_{inc.i} = \sqrt{I_{inc.i} + Q_{inc.i}}\\
        \varphi_{inc.i} = \arctan{(\frac{I_{inc.i}}{Q_{inc.i}})}
    \end{cases}
\end{equation}
$$
Дальше необходимо понять, что мы вычислили.


На основе модели известна необходимая разность фаз между приемными каналами, чтобы угол пеленга соответствовал истинном углу источника радиоизлучения (ИРИ). При решении задачи в первом приближении полагаем, что ошибки вносимые радиосканером несущественные и что ошибка определения углов связана исключительно с помеховым сигналом.

Надо понять от чего зависит дельта между должной разностью фаз и недолжной разностью фаз. Система становится уже сложной преобразуем систему уравнений (1)






