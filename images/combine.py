from PIL import Image

top = Image.open("margin_t_q16_n256_L64_Lc5_vf60.png")
bot = Image.open("../experiments/04_margin_tessellation/results/lkc_vs_delta_q16_n256.png")

scale = top.width / bot.width
new_h = int(bot.height * scale)
bot_r = bot.resize((top.width, new_h), Image.LANCZOS)

combined = Image.new("RGB", (top.width, top.height + new_h), (255, 255, 255))
combined.paste(top.convert("RGB"), (0, 0))
combined.paste(bot_r.convert("RGB"), (0, top.height))
combined.save("margin_t_q16_n256_L64_Lc5_vf60_combined.png")
print(f"Saved: {top.width}x{top.height + new_h}")
