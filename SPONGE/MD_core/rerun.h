#pragma once

#include <cstdint>

struct RERUN_information
{
    CONTROLLER* controller = NULL;
    MD_INFORMATION* md_info =
        NULL;  // 指向自己主结构体的指针，以方便调用主结构体的信息
    FILE* traj_file = NULL;
    FILE* box_file = NULL;
    FILE* vel_file = NULL;
    VECTOR* coordinate_staging = NULL;
    VECTOR* velocity_staging = NULL;
    bool has_velocity_stream = false;
    bool has_frame = false;
    VECTOR frame_box_length;
    VECTOR frame_box_angle;
    std::size_t box_line_number = 0;
    std::uint64_t complete_frames_read = 0;
    LTMatrix3 g;  // 盒子变化的速度
    int need_box_update = 0;
    int start_frame = 0;
    int strip_frame = 0;
    void Initial(CONTROLLER* controller, MD_INFORMATION* md_info);
    int Read_Frame(VECTOR* box_length, VECTOR* box_angle);
    void Close_Files();
    void Clear();
    bool Iteration(int strip = -1);
};
