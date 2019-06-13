#  Optimization direction

Through the previous introduction, we can find that there are two main directions to improve the effect. One is to improve the algorithm, and improve the speed of the program by improving the algorithms such as the rays generation algorithm, the hit algorithm, and the sampling algorithm. The other is to optimize the computational structure, such as parallelizing the program to speed up the operation.

These two directions are both normal solutions for real time rendering.But, we also find some different ways which have really good effect. 

NVIDIA published their real time rendering technology in 2018. A 24-frame 1080p video shocked everyone, which means that real-time rendering in the game is no longer out of reach. We are also curious how it is implemented.

Actually, it uses a new kind of Turing GPU with RT cores. Every RT core can provide 10 gaga rays per second processing speed and basic ray triangle intersection model. But  for us, more important is the software. It uses three structures to achieve more efficient ray tracing include RayTracing Pipeline,   Acceleration Structures and Shader Table . The idea is to decouple as much as possible between modules to improve parallelism. The RT core can provide fast traversal and intersection, and the shader table can achieve more realistic coloring effects. 

Unlike other methods, it reduces the sampling rate of ray tracing, and optimizes the resulting image using a denoiser and achieves good results. The shocking thing is that in the case of 1 sample per pixel, the picture after noise reduction is comparable to the real shot.

In summary, there are two ways to improve the ray tracing effect in the future. One is the conventional approach, improving the algorithm or increasing parallelism. The second is to first render and then denoise through other ways just as NVIDIA.

