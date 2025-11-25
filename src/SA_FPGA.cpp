#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <math.h>
#include <random>
#include <chrono>
#include <unordered_set>

std::mt19937 rng(123456);

using namespace std;
using ns = std::chrono::nanoseconds;
using Clock = std::chrono::high_resolution_clock;

int R = 0, C = 0;

struct pin_info
{
    string name;
    double x;
    double y;
};

struct logic_info
{
    string name;
    int x;
    int y;
    vector<int> connect_nets;
};
struct net_info
{
    string name;
    int degree;
    vector<string> connections;
};
struct rec
{
    int xmax;
    int xmin;
    int ymax;
    int ymin;
};
struct net_hpwl_cong
{
    int hpwl;
    rec box;
};
struct change_info
{
    string nameA, nameB;
    int Ax, Ay;
    int Bx, By;
};

void initial_place(vector<vector<string>> &array2D, int R, int C, vector<logic_info> &logic_blocks_table)
{
    int idx = 0;
    for (int i = 0; i < R && idx < (int)logic_blocks_table.size(); i++)
    {
        for (int j = 0; j < C && idx < (int)logic_blocks_table.size(); j++)
        {
            array2D[i][j] = logic_blocks_table[idx].name;
            logic_blocks_table[idx].x = i;
            logic_blocks_table[idx].y = j;
            idx++;
        }
    }
}

rec calculate_bounding_box(const net_info &net,
                           const vector<logic_info> &logic_blocks_table,
                           const vector<pin_info> &pin_table,
                           const unordered_map<string, int> &name2idx,
                           const unordered_map<string, int> &pin2idx)
{
    bool first = true;
    double xmin = 0.0, xmax = 0.0, ymin = 0.0, ymax = 0.0;

    for (const auto &item : net.connections)
    {
        auto it_blk = name2idx.find(item);
        auto it_pin = pin2idx.find(item);

        if (it_blk != name2idx.end())
        {
            // movable block
            int idx = it_blk->second;
            double xb = (double)logic_blocks_table[idx].x;
            double yb = (double)logic_blocks_table[idx].y;
            double this_xmin = xb;
            double this_xmax = xb + 1.0;
            double this_ymin = yb;
            double this_ymax = yb + 1.0;

            if (first)
            {
                xmin = this_xmin;
                xmax = this_xmax;
                ymin = this_ymin;
                ymax = this_ymax;
                first = false;
            }
            else
            {
                xmin = min(xmin, this_xmin);
                xmax = max(xmax, this_xmax);
                ymin = min(ymin, this_ymin);
                ymax = max(ymax, this_ymax);
            }
        }
        else if (it_pin != pin2idx.end())
        {
            // fixed I/O pin
            int idx = it_pin->second;
            double xp = pin_table[idx].x;
            double yp = pin_table[idx].y;
            double this_xmin = xp;
            double this_xmax = xp; // pin 是點
            double this_ymin = yp;
            double this_ymax = yp;

            if (first)
            {
                xmin = this_xmin;
                xmax = this_xmax;
                ymin = this_ymin;
                ymax = this_ymax;
                first = false;
            }
            else
            {
                xmin = min(xmin, this_xmin);
                xmax = max(xmax, this_xmax);
                ymin = min(ymin, this_ymin);
                ymax = max(ymax, this_ymax);
            }
        }
        // 若 item 既不是 block 又不是 pin，就忽略
    }

    // 萬一這個 net 沒有任何合法 terminal（理論上不會）
    if (first)
    {
        rec deg;
        deg.xmin = deg.xmax = deg.ymin = deg.ymax = 0.0;
        return deg;
    }

    rec Rect;
    Rect.xmin = xmin;
    Rect.xmax = xmax;
    Rect.ymin = ymin;
    Rect.ymax = ymax;
    return Rect;
}

void compute_hpwl_cong(const vector<logic_info> &logic_blocks_table,
                       const vector<pin_info> &pin_table,
                       const vector<net_info> &net_table,
                       const unordered_map<string, int> &name2idx,
                       const unordered_map<string, int> &pin2idx,
                       vector<vector<int>> &U,
                       unordered_map<string, net_hpwl_cong> &net2hpwl_cong,
                       double &total_hpwl,
                       double &CC)
{
    // 初始化 coverage
    for (int y = 0; y < R; ++y)
    {
        fill(U[y].begin(), U[y].end(), 0);
    }

    total_hpwl = 0.0;

    for (size_t i = 0; i < net_table.size(); ++i)
    {
        const net_info &net = net_table[i];
        rec box = calculate_bounding_box(net, logic_blocks_table, pin_table,
                                         name2idx, pin2idx);
        double hpwl = (box.xmax - box.xmin) + (box.ymax - box.ymin);
        total_hpwl += hpwl;

        net2hpwl_cong[net.name].hpwl = hpwl;
        net2hpwl_cong[net.name].box = box;
        int x_start = max(0, (int)floor(box.xmin));
        int x_end = min(C, (int)ceil(box.xmax)); // [x_start, x_end)
        int y_start = max(0, (int)floor(box.ymin));
        int y_end = min(R, (int)ceil(box.ymax)); // [y_start, y_end)

        for (int y = y_start; y < y_end; ++y)
        {
            for (int x = x_start; x < x_end; ++x)
            {
                U[y][x] += 1;
            }
        }
    }

    // 計算 CC
    long long N = (long long)R * (long long)C;
    long long sumU = 0;
    long long sumU2 = 0;
    for (int y = 0; y < R; ++y)
    {
        for (int x = 0; x < C; ++x)
        {
            int v = U[y][x];
            sumU += v;
            sumU2 += 1LL * v * v;
        }
    }

    if (sumU == 0)
    {
        CC = 1.0; // 沒有任何 net 覆蓋時，視為完美均勻
    }
    else
    {
        CC = (double)N * (double)sumU2 / ((double)sumU * (double)sumU);
    }
}

double compute_cost(double total_hpwl, double CC, double lambda)
{
    /*const double CC_target = 1.12;              // 你希望的 CC 水準
    double pen = std::max(0.0, CC - CC_target); // 只罰超過的那一段*/
    return total_hpwl * CC;
}

change_info random_move(vector<logic_info> &logic_blocks_table, unordered_map<string, int> &name2idx,
                        vector<vector<string>> &array2D,
                        int R, int C) // 原來 new_x,new_y 放什麼（block name 或 "."）
{
    int n = logic_blocks_table.size();
    std::uniform_int_distribution<int> dist_block(0, n - 1);
    std::uniform_int_distribution<int> distR(0, R - 1);
    std::uniform_int_distribution<int> distC(0, C - 1);
    change_info changes;

    // 隨機選一個 block
    int idx = dist_block(rng);

    auto &A = logic_blocks_table[idx];
    int Ax = logic_blocks_table[idx].x;
    int Ay = logic_blocks_table[idx].y;

    // 隨機選一個 grid cell（可能是空格、可能是 block）
    int Bx = distR(rng);
    int By = distC(rng);

    // 若選到原地（不動），就重選
    while (Bx == Ax && By == Ay)
    {
        Bx = distR(rng);
        By = distC(rng);
    }

    // 記錄要被換走的（可能是 block name，也可能是 "."）
    string swap_target = array2D[Bx][By];
    string temp;
    // ========== Case 1: target 是 block → swap ==========
    if (swap_target != "none")
    {
        int idxB = name2idx[swap_target];
        auto &B = logic_blocks_table[idxB];
        temp = B.name;
        A.x = Bx;
        A.y = By;
        B.x = Ax;
        B.y = Ay;

        // 更新 array2D
        array2D[Ax][Ay] = swap_target;
        array2D[Bx][By] = A.name;
        /*cout << "[RANDOM MOVE] "
         << "A = " << A.name << " from (" << Ax << ", " << Ay << ")"
         << "  <-->  B = " << swap_target
         << " at (" << Bx << ", " << By << ")\n";*/
    }
    // ========== Case 2: target 是空格 → 單獨移動 ==========
    else
    {
        A.x = Bx;
        A.y = By;
        temp = "none";
        array2D[Bx][By] = A.name;
        array2D[Ax][Ay] = "none";
        /*cout << "[RANDOM MOVE] "
         << "A = " << A.name << " from (" << Ax << ", " << Ay << ")"
         << "  -->  EMPTY at (" << Bx << ", " << By << ")\n";*/
    }
    changes.nameA = A.name;
    changes.Ax = Ax;
    changes.Ay = Ay;
    changes.nameB = temp;
    changes.Bx = Bx;
    changes.By = By;
    return changes;
}

void undo_move(vector<logic_info> &logic_blocks_table,
               unordered_map<string, int> &name2idx,
               vector<vector<string>> &array2D,
               const change_info &chg)
{
    // 還原 A
    int idxA = name2idx[chg.nameA];
    auto &A = logic_blocks_table[idxA];

    if (chg.nameB != "none")
    {
        // block <-> block 的情況
        int idxB = name2idx[chg.nameB];
        auto &B = logic_blocks_table[idxB];

        A.x = chg.Ax;
        A.y = chg.Ay;
        B.x = chg.Bx;
        B.y = chg.By;

        array2D[chg.Ax][chg.Ay] = chg.nameA;
        array2D[chg.Bx][chg.By] = chg.nameB;
    }
    else
    {
        // block <-> 空格 的情況
        A.x = chg.Ax;
        A.y = chg.Ay;

        array2D[chg.Ax][chg.Ay] = chg.nameA;
        array2D[chg.Bx][chg.By] = "none";
    }
}

void SA(vector<logic_info> &logic_blocks_table,
        vector<vector<string>> &array2D,
        vector<pin_info> &pin_table,
        vector<net_info> &net_table,
        unordered_map<string, int> &name2idx,
        unordered_map<string, int> &pin2idx,
        unordered_map<string, net_hpwl_cong> &net2hpwl_cong,
        vector<vector<int>> &U,
        double &total_hpwl,
        double &CC)
{
    double pen0 = std::max(0.0, CC - 1.0);
    if (pen0 < 1e-6)
        pen0 = 1.0; // 避免除 0
    double lambda = 0.15 * total_hpwl / pen0;
    double cur_cost = compute_cost(total_hpwl, CC, lambda);
    double best_cost = cur_cost;

    cout << "lambda: " << lambda << endl;

    // 紀錄最佳解
    auto best_logic_blocks = logic_blocks_table;
    auto best_array2D = array2D;
    double best_hpwl = total_hpwl;
    double best_CC = CC;

    cerr << "[INIT] HPWL=" << total_hpwl
         << " CC=" << CC
         << " cost=" << cur_cost << endl;

    // SA 參數
    double T = 50.0;
    double T_end = 1e-5;
    double alpha = 0.95;
    int iter_per_T = 5000;

    std::uniform_real_distribution<double> dist01(0.0, 1.0);

    // 先從現在的 U 算一次 sumU / sumU2，之後都 incremental 更新
    const long long N = (long long)R * (long long)C;
    long long sumU = 0;
    long long sumU2 = 0;
    for (int y = 0; y < R; ++y)
    {
        for (int x = 0; x < C; ++x)
        {
            int v = U[y][x];
            sumU += v;
            sumU2 += 1LL * v * v;
        }
    }

    // 給 incremental 用的小結構
    struct NetUpdate
    {
        int net_idx;
        std::string net_name;
        double old_hpwl;
        rec old_box;
        double new_hpwl;
        rec new_box;
    };

    struct CellBackup
    {
        int x, y;
        int old_val;
    };

    while (T > T_end)
    {
        for (int it = 0; it < iter_per_T; ++it)
        {
            // 1. 做一個 random move（可能 block<->block 或 block<->空格）
            change_info chg = random_move(logic_blocks_table, name2idx,
                                          array2D, R, C);

            // 2. 找出受影響的 nets（A、B 兩顆 cell 相關聯的 nets）
            std::unordered_set<int> affected_net_set;

            auto collect_nets = [&](const std::string &blk_name)
            {
                if (blk_name == "none")
                    return;
                auto it_blk = name2idx.find(blk_name);
                if (it_blk == name2idx.end())
                    return;
                const auto &blk = logic_blocks_table[it_blk->second];
                for (int nid : blk.connect_nets)
                {
                    affected_net_set.insert(nid);
                }
            };

            collect_nets(chg.nameA);
            collect_nets(chg.nameB);

            // 3. 只對 affected nets 用 net2hpwl_cong 查舊 hpwl / box，
            //    然後重新算新的 bbox / hpwl
            std::vector<NetUpdate> updates;
            updates.reserve(affected_net_set.size());

            double trial_hpwl = total_hpwl; // 從現在的 total_hpwl 開始加差值

            for (int nid : affected_net_set)
            {
                const net_info &net = net_table[nid];
                const std::string &net_name = net.name;

                // 從 net2hpwl_cong 查舊 hpwl / box
                auto it_info = net2hpwl_cong.find(net_name);
                double old_hp = 0.0;
                rec old_rect{};
                if (it_info != net2hpwl_cong.end())
                {
                    old_hp = (double)it_info->second.hpwl;
                    old_rect = it_info->second.box;
                }

                // 用目前 (已 move) 的 logic_blocks_table/pin_table/name2idx/pin2idx
                // 重新算 bbox
                rec new_rect = calculate_bounding_box(net,
                                                      logic_blocks_table,
                                                      pin_table,
                                                      name2idx,
                                                      pin2idx);
                double new_hp = (new_rect.xmax - new_rect.xmin) +
                                (new_rect.ymax - new_rect.ymin);

                trial_hpwl += (new_hp - old_hp);

                NetUpdate nu;
                nu.net_idx = nid;
                nu.net_name = net_name;
                nu.old_hpwl = old_hp;
                nu.old_box = old_rect;
                nu.new_hpwl = new_hp;
                nu.new_box = new_rect;
                updates.push_back(nu);
            }

            // 4. 對 U 做「暫時性的」更新，順便算 new_sumU / new_sumU2
            long long new_sumU = sumU;
            long long new_sumU2 = sumU2;

            std::vector<CellBackup> cell_backups;
            cell_backups.reserve(256);
            std::unordered_map<long long, int> visited; // key: (y<<32)|x -> index in backups

            auto apply_rect_delta = [&](const rec &box, int delta)
            {
                int x_start = std::max(0, (int)floor((double)box.xmin));
                int x_end = std::min(C, (int)ceil((double)box.xmax)); // [x_start, x_end)
                int y_start = std::max(0, (int)floor((double)box.ymin));
                int y_end = std::min(R, (int)ceil((double)box.ymax)); // [y_start, y_end)

                for (int y = y_start; y < y_end; ++y)
                {
                    for (int x = x_start; x < x_end; ++x)
                    {
                        long long key = ((long long)y << 32) | (unsigned int)x;
                        int idx_bk;
                        auto itv = visited.find(key);
                        if (itv == visited.end())
                        {
                            CellBackup bk{x, y, U[y][x]};
                            cell_backups.push_back(bk);
                            idx_bk = (int)cell_backups.size() - 1;
                            visited[key] = idx_bk;
                        }
                        else
                        {
                            idx_bk = itv->second;
                        }

                        int old_v = U[y][x];
                        int new_v = old_v + delta;
                        U[y][x] = new_v;

                        new_sumU += delta;
                        new_sumU2 += 1LL * new_v * new_v - 1LL * old_v * old_v;
                    }
                }
            };

            // 用 net2hpwl_cong 的舊 box 先 -1，再用新 box +1
            for (const auto &u : updates)
            {
                if (u.old_hpwl > 0.0)
                {
                    apply_rect_delta(u.old_box, -1);
                }
                apply_rect_delta(u.new_box, +1);
            }

            // 用更新後的 U 統計值得到 new_CC
            double new_CC;
            if (new_sumU == 0)
            {
                new_CC = 1.0;
            }
            else
            {
                new_CC = (double)N * (double)new_sumU2 /
                         ((double)new_sumU * (double)new_sumU);
            }

            double new_cost = compute_cost(trial_hpwl, new_CC, lambda);
            double delta_cost = new_cost - cur_cost;
            bool accept = false;

            // 5. SA 接受規則
            if (delta_cost <= 0)
            {
                accept = true;
            }
            else
            {
                double u = dist01(rng);
                if (u < exp(-delta_cost / T))
                    accept = true;
            }

            if (accept)
            {
                // 6. 接受：commit incremental 更新
                cur_cost = new_cost;
                total_hpwl = trial_hpwl;
                CC = new_CC;
                sumU = new_sumU;
                sumU2 = new_sumU2;

                // 用 net2hpwl_cong 更新 affected nets 的 hpwl / box
                for (const auto &u : updates)
                {
                    auto &info = net2hpwl_cong[u.net_name];
                    info.hpwl = (int)u.new_hpwl;
                    info.box = u.new_box;
                }

                // 更新 best solution
                if (cur_cost < best_cost)
                {
                    best_cost = cur_cost;
                    best_hpwl = total_hpwl;
                    best_CC = CC;
                    best_logic_blocks = logic_blocks_table;
                    best_array2D = array2D;
                }
            }
            else
            {
                // 7. 不接受：把 U 還原，位置用 undo_move 拉回
                for (const auto &bk : cell_backups)
                {
                    U[bk.y][bk.x] = bk.old_val;
                }
                undo_move(logic_blocks_table, name2idx, array2D, chg);
                // total_hpwl / CC / sumU / sumU2 / net2hpwl_cong 都維持舊的
            }
        }

        cerr << "[TEMP] T=" << T
             << " best_cost=" << best_cost
             << " best_HPWL=" << best_hpwl
             << " best_CC=" << best_CC << endl;

        T *= alpha; // 降溫
    }

    // 退火結束 → 回到最佳解
    logic_blocks_table = best_logic_blocks;
    array2D = best_array2D;
    total_hpwl = best_hpwl;
    CC = best_CC;

    // 最後再 full recompute 一次，讓 net2hpwl_cong / U 對應到 best solution
    compute_hpwl_cong(logic_blocks_table, pin_table, net_table,
                      name2idx, pin2idx,
                      U, net2hpwl_cong, total_hpwl, CC);
    best_cost = compute_cost(total_hpwl, CC, lambda);

    cerr << "[FINAL] HPWL=" << total_hpwl
         << " CC=" << CC
         << " cost=" << best_cost << endl;
}

void parse_input(ifstream &parse, vector<logic_info> &logic_blocks_table, vector<pin_info> &pin_table,
                 vector<net_info> &net_table, unordered_map<string, int> &name2idx, unordered_map<string, int> &pin2idx, unordered_map<string, int> &net2idx)
{
    string line;
    getline(parse, line);
    istringstream info(line);
    int logic_block_R, logic_block_C, num_logic_blocks, num_io_pins, num_nets;
    info >> logic_block_R >> logic_block_C >> num_logic_blocks >> num_io_pins >> num_nets;
    R = logic_block_R;
    C = logic_block_C;
    for (int i = 0; i < num_logic_blocks; i++)
    {
        string name;
        getline(parse, line);
        istringstream logic_inf(line);
        logic_inf >> name;
        logic_info inst;
        inst.name = name;
        inst.x = 0;
        inst.y = 0;
        logic_blocks_table.push_back(inst);
        name2idx.insert({name, i});
    }
    for (int i = 0; i < num_io_pins; i++)
    {
        string name;
        double x, y;
        getline(parse, line);
        istringstream pin_inf(line);
        pin_inf >> name >> x >> y;
        pin_info inst;
        inst.name = name;
        inst.x = x;
        inst.y = y;
        pin_table.push_back(inst);
        pin2idx.insert({name, i});
    }
    for (int i = 0; i < num_nets; i++)
    {
        string name;
        int degree;
        getline(parse, line);
        istringstream net_inf(line);
        net_inf >> name >> degree;
        net_info inst;
        inst.name = name;
        inst.degree = degree;
        for (int j = 0; j < degree; j++)
        {
            string item;
            net_inf >> item;
            inst.connections.push_back(item);
            auto it = name2idx.find(item);
            if (it != name2idx.end())
            {
                logic_blocks_table[it->second].connect_nets.push_back(i);
            }
        }
        net_table.push_back(inst);
        net2idx.insert({name, i});
    }
}

void print_logic_blocks(const vector<logic_info> &logic_blocks_table)
{
    cout << "=== Logic Blocks ===\n";
    for (size_t i = 0; i < logic_blocks_table.size(); ++i)
    {
        const auto &lb = logic_blocks_table[i];
        cout << "[" << i << "] "
             << "name = " << lb.name
             << ", pos = (" << lb.x << ", " << lb.y << ")"
             << ", nets = { ";

        for (size_t k = 0; k < lb.connect_nets.size(); ++k)
        {
            cout << lb.connect_nets[k];
            if (k + 1 != lb.connect_nets.size())
                cout << ", ";
        }

        cout << " }\n";
    }
    cout << endl;
}

void print_pins(const vector<pin_info> &pin_table)
{
    cout << "=== IO Pins ===\n";
    for (size_t i = 0; i < pin_table.size(); ++i)
    {
        const auto &p = pin_table[i];
        cout << i << ": name=" << p.name
             << " x=" << p.x
             << " y=" << p.y << "\n";
    }
    cout << endl;
}
void print_nets(const vector<net_info> &net_table)
{
    cout << "=== Nets ===\n";
    for (size_t i = 0; i < net_table.size(); ++i)
    {
        const auto &n = net_table[i];
        cout << i << ": name=" << n.name
             << " degree=" << n.degree
             << " connections = ";
        for (auto &c : n.connections)
        {
            cout << c << " ";
        }
        cout << "\n";
    }
    cout << endl;
}

void print_net2hpwl_cong(const unordered_map<string, net_hpwl_cong> &net2hpwl_cong)
{
    cout << "=== net2hpwl_cong ===\n";
    int count = 0;
    for (const auto &kv : net2hpwl_cong)
    {
        const string &net_name = kv.first;
        const net_hpwl_cong &info = kv.second;

        cout << "Net: " << net_name << "\n";
        cout << "  HPWL = " << info.hpwl << "\n";
        cout << "  Box: "
             << "(xmin=" << info.box.xmin
             << ", xmax=" << info.box.xmax
             << ", ymin=" << info.box.ymin
             << ", ymax=" << info.box.ymax
             << ")\n";
        count++;
    }
    cout << count << endl;
    cout << endl;
}

void print_U(const vector<vector<int>> &U)
{
    cout << "=== U congestion map ===\n";
    for (int r = 0; r < (int)U.size(); r++)
    {
        for (int c = 0; c < (int)U[r].size(); c++)
        {
            cout << U[r][c] << " ";
        }
        cout << "\n";
    }
    cout << endl;
}

void print_array2D(const vector<vector<string>> &array2D)
{
    cout << "=== array2D placement ===\n";
    for (int r = 0; r < (int)array2D.size(); r++)
    {
        for (int c = 0; c < (int)array2D[r].size(); c++)
        {
            cout << array2D[r][c] << " ";
        }
        cout << "\n";
    }
    cout << endl;
}

void write_output(const string &out_file,
                  const vector<logic_info> &logic_blocks_table)
{
    ofstream fout(out_file);
    if (!fout)
    {
        cerr << "Failed to open output file: " << out_file << endl;
        return;
    }
    for (const auto &lb : logic_blocks_table)
    {
        fout << lb.name << " " << lb.x << " " << lb.y << "\n";
    }
    fout.close();
}

int main(int argc, char *argv[])
{
    auto start = Clock::now();

    if (argc < 3)
    {
        cerr << "Usage: " << argv[0] << " <lef_file> <def_file>" << endl;
        return 1;
    }

    ifstream parse_in(argv[1]);
    if (!parse_in.is_open())
    {
        cerr << "Failed to open input file: " << argv[1] << endl;
        return 1;
    }
    vector<logic_info> logic_blocks_table;
    vector<pin_info> pin_table;
    vector<net_info> net_table;
    unordered_map<string, int> name2idx;
    unordered_map<string, int> pin2idx;
    unordered_map<string, int> net2idx;
    unordered_map<string, net_hpwl_cong> net2hpwl_cong;
    double total_hpwl, CC;
    parse_input(parse_in, logic_blocks_table, pin_table, net_table, name2idx, pin2idx, net2idx);
    vector<vector<string>> array2D(R, vector<string>(C, "none"));
    vector<vector<int>> U(R, vector<int>(C, 0));
    initial_place(array2D, R, C, logic_blocks_table);
    // print_array2D(array2D);
    //  Optional: print parsed data for verification
    compute_hpwl_cong(logic_blocks_table, pin_table, net_table, name2idx, pin2idx, U, net2hpwl_cong, total_hpwl, CC);
    SA(logic_blocks_table, array2D,
       pin_table, net_table,
       name2idx, pin2idx,
       net2hpwl_cong, U,
       total_hpwl, CC);

    // 看一下最後的 U / net2hpwl_cong（可選）
    // print_net2hpwl_cong(net2hpwl_cong);
    // print_U(U);
    // print_array2D(array2D);

    cout << "[FINAL] HPWL=" << total_hpwl
         << " CC=" << CC << endl;
    write_output(string(argv[2]), logic_blocks_table);

    // print_net2hpwl_cong(net2hpwl_cong);
    // print_U(U);
    // cout << total_hpwl << " " << CC << endl;
    //  print_logic_blocks(logic_blocks_table);
    //  print_pins(pin_table);
    //  print_nets(net_table);
    auto end = Clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    cout << "Execution time: " << duration.count() << " milliseconds" << endl;

    return 0;
}