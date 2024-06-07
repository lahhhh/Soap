#pragma once

#include <concepts>
#include <type_traits>
#include <ranges>

#include <QMap>

namespace custom {

	template <typename T>
	concept is_integral = std::is_integral_v<T>;

	template <typename T1, typename T2>
	concept not_same = !std::same_as<T1, T2>;

	/*     element type       */

	template <typename T>
	using get_element_raw_type = std::decay_t<decltype(*std::begin(std::declval<std::decay_t<T>>()))>;

	template <typename T>
	using get_element_arg_type = decltype(*std::begin(std::declval<T>()));

	template <typename T1, typename T2>
	struct change_element_type_helper;

	template <template<typename...>typename Outer, typename Original, typename To, typename... Ts>
	struct change_element_type_helper<Outer<Original, Ts...>, To> {

		using type = Outer<To>;
	};

	template <template<typename, int, int...>typename Outer, typename Original, typename To, int Int, int... Ints>
	struct change_element_type_helper<Outer<Original, Int, Ints...>, To> {

		using type = Outer<To, Int, Ints...>;
	};

	template <typename Key, typename Value, typename To>
	struct change_element_type_helper<QMap<Key, Value>, To> {

		using type = QMap<Key, To>;
	};

	template <typename T, typename ToEle>
	using change_element_type = change_element_type_helper<T, ToEle>::type;

	template <template<typename...>typename Outer, typename T>
	struct change_outer_type_helper {
		using Ele = get_element_raw_type<T>;
		using type = Outer<Ele>;
	};

	/*     return type deduction     */
	template <typename Fun, typename ...Args>
	using get_return_type = decltype(std::declval<Fun>()(std::declval<Args...>()));

	template <typename Fun, typename ...Args>
	struct sapply_return_type_helper {
		using return_type = QList<get_return_type<Fun, Args...>>;
	};

	template <typename Fun, typename ...Args>
		requires std::same_as<get_return_type<Fun, Args...>, void>
	struct sapply_return_type_helper<Fun, Args...> {
		using return_type = void;
	};

	template <typename Fun, typename ...Args>
	using sapply_return_type = sapply_return_type_helper<Fun, Args...>::return_type;

	template <typename T>
	concept less_comparable = requires(const std::remove_reference_t<T>&x, const std::remove_reference_t<T>&y) { {x < y} -> std::same_as<bool>; };

	template <typename T1, typename T2>
	concept same_original_type = std::same_as<std::decay_t<T1>, std::decay_t<T2>>;

	template <typename T>
	concept is_container = requires (const std::decay_t<T>&container) {
		{ std::cbegin(container) } -> std::forward_iterator;
		{ std::cend(container) } /*-> std::sentinel_for<decltype(std::cbegin(container))>*/;
		{ std::ranges::size(container) } -> std::convertible_to<std::size_t>;
	};

	template <typename T>
	concept not_container = !is_container<T>;

	template <typename T, typename Ele>
	concept is_specific_container = is_container<T> && std::same_as<get_element_raw_type<T>, Ele>;

	template <typename T>
	concept is_slice_container = is_specific_container<T, bool>;

	template <typename T>
	concept is_bool_container = is_slice_container<T>;

	template <typename Order>
	concept is_order_container = is_container<Order>
		&& std::is_integral_v<get_element_raw_type<Order>>
		&& !std::same_as<get_element_raw_type<Order>, bool>;

	template <typename T>
	concept is_arithmetic_container = is_container<T> && std::is_arithmetic_v<get_element_raw_type<T> >;

	template <typename T>
	concept random_accessible = requires (const std::decay_t<T>& container) {
		{container[std::declval<const std::size_t&>()]};
	};

	template <typename T>
	concept is_random_access_container = is_container<T> && random_accessible<T>;

	template <typename T>
	concept construct_from_size = std::is_constructible_v<T, std::size_t>;

	template <typename T1, typename T2>
	concept same_element = is_container<T1> && is_container<T2> && std::same_as<get_element_raw_type<T1>, get_element_raw_type<T2>>;
}